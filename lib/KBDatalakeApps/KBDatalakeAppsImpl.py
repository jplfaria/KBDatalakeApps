# -*- coding: utf-8 -*-
#BEGIN_HEADER
import json
import logging
import os
import uuid
from pathlib import Path
import shutil
import subprocess
import time
import polars as pl
import cobra
import pandas as pd

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.RAST_SDKClient import RAST_SDK
from installed_clients.kb_baktaClient import kb_bakta
from installed_clients.kb_psortbClient import kb_psortb
from installed_clients.kb_kofamClient import kb_kofam
from modelseedpy import MSModelUtil
from modelseedpy.core.msgenome import MSGenome, MSFeature
from cobrakbase import KBaseAPI
from installed_clients.baseclient import ServerError
from annotation.annotation import test_annotation, run_rast, run_kofam
from executor.task_executor import TaskExecutor
from executor.task import task_rast, task_kofam, task_psortb
from KBDatalakeApps.KBDatalakeUtils import KBDataLakeUtils,run_phenotype_simulation,run_model_reconstruction

# Import KBUtilLib utilities for common functionality
#from kbutillib import KBWSUtils, KBCallbackUtils, SharedEnvUtils

#class DatalakeAppUtils(KBWSUtils, KBCallbackUtils, SharedEnvUtils):
#    """Custom utility class combining KBUtilLib modules for datalake operations."""
#    pass

def human_size(size):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024:
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{size:.1f}PB"


def print_path(root: Path):
    if not root.exists():
        print(f"{root} does not exist")
        return

    print(root.name)
    _print_tree(root, prefix="")


def _print_tree(root: Path, prefix: str):
    # skip folder mmseqs2_tmp
    if root.name == 'mmseqs2_tmp':
        return
    entries = sorted(root.iterdir(), key=lambda p: (p.is_file(), p.name.lower()))
    for i, path in enumerate(entries):
        is_last = i == len(entries) - 1
        connector = "└── " if is_last else "├── "

        if path.is_file():
            size = human_size(path.stat().st_size)
            print(prefix + connector + f"{path.name} ({size})")
        else:
            print(prefix + connector + path.name)
            extension = "    " if is_last else "│   "
            _print_tree(path, prefix + extension)


def get_berdl_token():
    return os.environ.get('KBASE_SECURE_CONFIG_PARAM_kbaselakehouseserviceaccount_token')
#END_HEADER


class KBDatalakeApps:
    '''
    Module Name:
    KBDatalakeApps

    Module Description:
    A KBase module: KBDatalakeApps

This module provides applications for interacting with the KBase data lake,
including data retrieval, processing, and analysis utilities.

Author: chenry
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:kbaseapps/KBDatalakeApps.git"
    GIT_COMMIT_HASH = "57ad9974af809ff24506351db02203a0481358e4"

    #BEGIN_CLASS_HEADER
    def _validate_params(self, params, required_keys):
        """Validate that required parameters are present."""
        for key in required_keys:
            if key not in params or params[key] is None:
                raise ValueError(f"Required parameter '{key}' is missing")

    @staticmethod
    def run_genome_pipeline(input_file):
        cmd = ["/kb/module/scripts/run_genome_pipeline.sh", str(input_file)]

        env = os.environ.copy()
        env.pop("PYTHONPATH", None)

        process = subprocess.Popen(
            cmd,
            stdout=None,  # inherit parent stdout
            stderr=None,  # inherit parent stderr
            env=env
        )

        ret = process.wait()
        if ret != 0:
            raise RuntimeError(
                f"Genome pipeline failed with exit code {ret}"
            )

    @staticmethod
    def run_pangenome_pipeline(input_file, selected_member_id):
        cmd = ["/kb/module/scripts/run_pangenome_pipeline.sh", str(input_file), str(selected_member_id)]

        env = os.environ.copy()
        env.pop("PYTHONPATH", None)

        process = subprocess.Popen(
            cmd,
            stdout=None,  # inherit parent stdout
            stderr=None,  # inherit parent stderr
            env=env
        )

        ret = process.wait()
        if ret != 0:
            raise RuntimeError(
                f"Genome pipeline failed with exit code {ret}"
            )

    @staticmethod
    def run_RAST_annotation(input_filepath, output_filename, rast_client):
        sequence_hash = {}
        current_id = None
        current_seq = []
        with open(input_filepath) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        sequence_hash[current_id] = "".join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                elif line:
                    current_seq.append(line)
        if current_id:
            sequence_hash[current_id] = "".join(current_seq)
        proteins = []
        ids = []
        for id, sequence in sequence_hash.items():
            proteins.append(sequence)
            ids.append(id)

        annotate_protein_params = {'proteins': proteins}
        print('annotate_protein_params:', annotate_protein_params)

        result = rast_client.annotate_proteins(annotate_protein_params)
        print('rast annotation result', result)
        functions_list = result.get('functions', [])
        records = []
        for id, functions in zip(ids, functions_list):
            records.append({
                'id': id,
                'functions': functions
            })
        df = pd.DataFrame.from_records(records)
        df.to_csv(output_filename, sep='\t', index=False)
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.config = config
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.logger = logging.getLogger(__name__)

        # Initialize KBUtilib utilities
        self.dfu = DataFileUtil(self.callback_url)
        self.kbase_api = KBaseAPI(os.environ['KB_AUTH_TOKEN'], config=config)
        self.kb_bakta = kb_bakta(self.callback_url, service_ver='beta')
        self.kb_psortb = kb_psortb(self.callback_url, service_ver='beta')
        self.kb_kofam = kb_kofam(self.callback_url, service_ver='beta')
        self.rast_client = RAST_SDK(self.callback_url, service_ver='beta')
        self.util = KBDataLakeUtils(kbendpoint=config["kbase-endpoint"],reference_path="/data/",module_path="/kb/module")
        self.util.set_token(get_berdl_token(), namespace="berdl")
        print('polars thread pool', pl.thread_pool_size())
        #END_CONSTRUCTOR
        pass


    def build_genome_datalake_tables(self, ctx, params):
        """
        Build genome datalake tables from Genome or GenomeSet objects
        This function takes a list of Genome or GenomeSet references and builds
        datalake tables from them. Optionally saves generated models to the workspace.
        :param params: instance of type "BuildGenomeDatalakeTablesParams"
           (Parameters for building genome datalake tables input_refs - list
           of workspace references to Genome or GenomeSet objects suffix -
           string suffix to append to generated table names save_models -
           boolean (0/1) indicating if generated models should be saved
           workspace_name - the name of the workspace for saving results) ->
           structure: parameter "input_refs" of list of String, parameter
           "suffix" of String, parameter "save_models" of Long, parameter
           "workspace_name" of String
        :returns: instance of type "ReportResults" (Standard report output
           structure used by KBase apps) -> structure: parameter
           "report_name" of String, parameter "report_ref" of String,
           parameter "workspace" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN build_genome_datalake_tables
        self.logger.info(f"Building genome datalake tables with params: {params}")
        skip_annotation = params['skip_annotation'] == 1
        skip_pangenome = params['skip_pangenome'] == 1
        skip_genome_pipeline = params['skip_genome_pipeline'] == 1
        skip_modeling_pipeline = params['skip_modeling_pipeline'] == 1

        input_params = Path(self.shared_folder) / 'input_params.json'
        print(str(input_params.resolve()))
        with open(str(input_params.resolve()), 'w') as fh:
            _params = dict(params)
            _params['_ctx'] = ctx
            _params['_config'] = self.config

            print('to create a copy for debug:', _params)

            fh.write(json.dumps(_params))

        print(os.environ)
        #print('ctx', ctx)
        #print('contig', self.config)
        print('data dir')
        print(os.listdir('/data'))
        if os.path.exists('/data') and os.path.exists('/data/reference_data'):
            print(os.listdir('/data/reference_data'))

            if os.path.exists('/data/reference_data/berdl_db'):
                print('berdl:', os.listdir('/data/reference_data/berdl_db'))

        if not skip_annotation:
            test_annotation(self.kb_kofam, self.kb_bakta, self.kb_psortb, self.rast_client)

        #print('BERDL Token')
        #print(self.get_berdl_token())

        # Validate required parameters
        self._validate_params(params, ['input_refs', 'workspace_name'])

        workspace_name = params['workspace_name']
        input_refs = params['input_refs']
        suffix = params.get('suffix', ctx['token'])
        save_models = params.get('save_models', 0)

        if not skip_genome_pipeline:
            self.run_genome_pipeline(input_params.resolve())
        else:
            print('skip genome pipeline')

        executor = TaskExecutor(max_workers=4)
        path_user_genome = Path(self.shared_folder) / "genome"
        path_user_genome.mkdir(parents=True, exist_ok=True)
        tasks = []
        for filename_faa in os.listdir(str(path_user_genome)):
            if filename_faa.endswith('.faa'):
                print('found', filename_faa)
                if skip_annotation:
                    print('skip_annotation')
                else:
                    tasks.append(executor.run_task(task_rast,
                                                   path_user_genome / filename_faa,
                                                   self.rast_client))
                    tasks.append(executor.run_task(task_kofam,
                                                   path_user_genome / filename_faa,
                                                   self.kb_kofam))

        print('Task set barrier')
        for t in tasks:
            print(f'await for {t.args} {t.status}')
            t.wait()
        for t in tasks:
            print(t.status)
            print(t.result)
            print(t.traceback)

        path_pangenome = Path(self.shared_folder) / "pangenome"
        path_pangenome.mkdir(parents=True, exist_ok=True)
        for folder_pangenome in os.listdir(str(path_pangenome)):
            if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                print(f'Found pangenome folder: {folder_pangenome}')
                # run pangenome pipeline for - folder_pangenome
                if not skip_pangenome:
                    self.run_pangenome_pipeline(input_params.resolve(), folder_pangenome)
                else:
                    print('skip pangenome')

        tasks_pangeome = []
        path_pangenome = Path(self.shared_folder) / "pangenome"
        path_pangenome.mkdir(parents=True, exist_ok=True)
        for folder_pangenome in os.listdir(str(path_pangenome)):
            if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                print(f'Found pangenome folder: {folder_pangenome}')
                path_pangenome_members = path_pangenome / folder_pangenome / 'genome'
                if path_pangenome_members.exists():
                    for _f in os.listdir(str(path_pangenome_members)):
                        if _f.endswith('.faa'):
                            tasks_pangeome.append(executor.run_task(task_rast,
                                                                    path_pangenome_members / _f,
                                                                    self.rast_client))


                    """
                    try:
                        print(f"run kb_bakta annotation for {genome}")
                        self.logger.info(f"run annotation for {genome}")
                        start_time = time.perf_counter()
                        result = self.kb_bakta.annotate_proteins(proteins)
                        end_time = time.perf_counter()
                        print(f"Execution time: {end_time - start_time} seconds")
                        print(f'received results of type {type(result)} and size {len(result)}')
                    except Exception as ex:
                        print(f'nope {ex}')

                    try:
                        print(f"run kb_psortb annotation for {genome}")
                        self.logger.info(f"run annotation for {genome}")
                        start_time = time.perf_counter()
                        result = self.kb_psortb.annotate_proteins(proteins, "-n")
                        end_time = time.perf_counter()
                        print(f"Execution time: {end_time - start_time} seconds")
                        print(f'received results of type {type(result)} and size {len(result)}')
                    except Exception as ex:
                        print(f'nope {ex}')
                    """

        if not skip_modeling_pipeline:
            for input_ref in input_refs:
                info = self.util.get_object_info(input_ref)
                genome_tsv_path = path_user_genome / f'{info[1]}_genome.tsv'
                self.util.run_user_genome_to_tsv(input_ref, genome_tsv_path)

                # Print head of genome TSV for testing
                print(f"=== Head of genome TSV: {genome_tsv_path} ===")
                with open(genome_tsv_path, 'r') as f:
                    for i, line in enumerate(f):
                        if i >= 5:
                            break
                        print(line.rstrip())
                print("=" * 50)

                # Run model reconstruction
                model_output_path = path_user_genome / f'{info[1]}_model'
                classifier_dir = Path('/kb/module/data')
                run_model_reconstruction(str(genome_tsv_path), str(model_output_path), str(classifier_dir),kbversion=self.util.kb_version)

                # Print head of model output for testing
                model_data_file = str(model_output_path) + "_data.json"
                print(f"=== Head of model data: {model_data_file} ===")
                with open(model_data_file, 'r') as f:
                    content = f.read(2000)
                    print(content[:2000])
                print("=" * 50)

                # Run phenotype simulation
                phenotype_output_path = path_user_genome / f'{info[1]}_phenotypes.json'
                cobra_model_path = str(model_output_path) + "_cobra.json"
                run_phenotype_simulation(cobra_model_path, str(phenotype_output_path),max_phenotypes=5,kbversion=self.util.kb_version)

                # Print head of phenotype output for testing
                print(f"=== Head of phenotype results: {phenotype_output_path} ===")
                with open(phenotype_output_path, 'r') as f:
                    content = f.read(2000)
                    print(content[:2000])
                print("=" * 50)
        else:
            print('skip modeling pipeline')

        #t_end_time = time.perf_counter()
        #print(f"Total Execution time annotation: {t_end_time - t_start_time} seconds")



        # Create KBaseFBA.GenomeDataLakeTables

        output_object = {
            'name': 'no_name',
            'description': 'fake',
            'genomeset_ref': '77057/3/1',
            'pangenome_data': [{
                'pangenome_id': 'super_fake',
                'pangenome_taxonomy': 'alien',
                'user_genomes': [],
                'datalake_genomes': [],
                'sqllite_tables_handle_ref': 'KBH_248173'
                # Chris sqlite KBH_248173
                # random e. coli? assembly KBH_248118
            }],
        }

        # Done with all tasks
        print('Task barrier')
        for t in tasks_pangeome:
            print(f'await for {t.args} {t.status}')
            t.wait()
        for t in tasks_pangeome:
            print(t.status)
            print(t.result)
            print(t.traceback)

        executor.shutdown()
        print_path(Path(self.shared_folder).resolve())

        #saved_object_info = self.kbase_api.save_object('fake_output',
        #                           params['workspace_name'],
        #                                   'KBaseFBA.GenomeDataLakeTables',
        #                                   output_object,
        #                                   meta={}
        #)

        #print(saved_object_info)

        # Create report with results
        report_client = KBaseReport(self.callback_url)

        output_directory = os.path.join(self.shared_folder, str(uuid.uuid4()))
        shutil.copytree('/kb/module/data/html', output_directory)

        # Write app-config.json so the DataTables Viewer knows which object to display
        # Use the first input reference as the UPA for the viewer
        app_config = {
            "upa": "76990/Test2"
        }
        app_config_path = os.path.join(output_directory, 'app-config.json')
        with open(app_config_path, 'w') as f:
            json.dump(app_config, f, indent=4)
        self.logger.info(f"Wrote app-config.json with UPA: {app_config['upa']}")

        shock_id = self.dfu.file_to_shock({
            'file_path': output_directory,
            'pack': 'zip'
        })['shock_id']

        self.logger.info(f"HTML directory contents: {os.listdir(output_directory)}")
        self.logger.info(f"Shock ID: {shock_id}")


        html_report = [{
            'shock_id': shock_id,
            'name': 'index.html',
            'label': 'BERDL Tables',
            'description': 'BERDL Table Viewer'
        }]

        report_params = {
            'message': 'message_in_app hi!',
            'warnings': ['example warning'],
            'workspace_name': params['workspace_name'],
            'objects_created': [],
            'html_links': html_report,
            'direct_html_link_index': 0,
            # 'html_window_height': int(params['report_height']),
            'html_window_height': 800,
        }

        report_info = report_client.create_extended_report(report_params)

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'workspace': workspace_name
        }
        #END build_genome_datalake_tables

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method build_genome_datalake_tables return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
