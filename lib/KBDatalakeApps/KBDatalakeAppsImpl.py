# -*- coding: utf-8 -*-
#BEGIN_HEADER
import json
import logging
import os
import uuid
from pathlib import Path
import shutil
import subprocess
import polars as pl
import pandas as pd
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.RAST_SDKClient import RAST_SDK
from installed_clients.kb_baktaClient import kb_bakta
from installed_clients.kb_psortbClient import kb_psortb
from installed_clients.kb_kofamClient import kb_kofam
from cobrakbase import KBaseAPI
from cobrakbase.AbstractHandleClient import AbstractHandle as HandleService
from installed_clients.baseclient import ServerError
from annotation.annotation import test_annotation, run_rast, run_kofam
from executor.task_executor import TaskExecutor
from executor.task import task_rast, task_kofam, task_psortb, task_bakta
from KBDatalakeApps.KBDatalakeUtils import KBDataLakeUtils, generate_ontology_tables
from KBDatalakeApps.utils import upload_blob_file, print_path, get_classifier, read_rast_as_genome
from KBDatalakeApps.utils import input_refs_to_genome_refs, validate_genome_refs


DEBUG_MODE = False


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
                f"genome pipeline failed with exit code {ret}"
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
                f"pangenome pipeline failed with exit code {ret}"
            )

    @staticmethod
    def run_annotation_pipeline(filename_faa):
        cmd = ["/kb/module/scripts/run_annotation.sh", str(filename_faa)]

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
                f"annotation pipeline failed with exit code {ret}"
            )

    @staticmethod
    def run_build_table(input_file, selected_member_id):
        cmd = ["/kb/module/scripts/run_generate_table.sh", str(input_file), str(selected_member_id)]

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
                f"table builder pipeline failed with exit code {ret}"
            )

    @staticmethod
    def run_model_pipeline(input_file):
        cmd = ["/kb/module/scripts/run_model_pipeline.sh", str(input_file)]

        env = os.environ.copy()
        #env.pop("PYTHONPATH", None)

        process = subprocess.Popen(
            cmd,
            stdout=None,  # inherit parent stdout
            stderr=None,  # inherit parent stderr
            env=env
        )

        ret = process.wait()
        if ret != 0:
            raise RuntimeError(
                f"model pipeline failed with exit code {ret}"
            )
    
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
        self.gfu = GenomeFileUtil(self.callback_url)
        self.kbase_api = KBaseAPI(os.environ['KB_AUTH_TOKEN'], config=config)
        self.kb_bakta = kb_bakta(self.callback_url, service_ver='beta')
        self.kb_psortb = kb_psortb(self.callback_url, service_ver='beta')
        self.kb_kofam = kb_kofam(self.callback_url, service_ver='beta')
        self.rast_client = RAST_SDK(self.callback_url, service_ver='beta')
        self.util = None
        self.hs = None

        if DEBUG_MODE:
            print('WARNING DEBUG MODE')
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

        path_token = Path.home() / ".kbase"
        path_token.mkdir(exist_ok=True)
        path_token = Path.home() / ".kbase" / "token"
        with open(path_token, 'w') as fh:
            fh.write(ctx['token'])
        self.util = KBDataLakeUtils(kbendpoint=self.config["kbase-endpoint"], reference_path="/data/",
                                    module_path="/kb/module",token=ctx['token'],
                                    dfu_client=self.dfu,gfu_client=self.gfu, callback_url=self.callback_url)
        self.util.set_token(get_berdl_token(), namespace="berdl")
        self.hs = HandleService(self.config["handle-service-url"], token=ctx['token'])
        skip_save_genome_annotation = params['skip_save_genome_annotation'] == 1
        skip_annotation = params['skip_annotation'] == 1
        skip_pangenome = params['skip_pangenome'] == 1
        skip_genome_pipeline = params['skip_genome_pipeline'] == 1
        skip_modeling_pipeline = params['skip_modeling_pipeline'] == 1
        export_all_data = params['export_all_content'] == 1
        export_pangenome_data = params['export_pangenome_data'] == 1
        export_genome_data = params['export_genome_data'] == 1
        export_databases = params['export_databases'] == 1
        export_folder_models = params['export_folder_models'] == 1
        export_folder_phenotypes = params['export_folder_phenotypes'] == 1
        input_refs = params['input_refs']
        output_object_name = params['output']
        suffix = params.get('suffix', ctx['token'])  # FIXME: why ctx token ???
        save_models = params.get('save_models', 0)
        workspace_name = params['workspace_name']
        param_clade_limit = 5  # max number of clades to process

        input_params = Path(self.shared_folder) / 'input_params.json'
        print(str(input_params.resolve()))

        genome_refs_all = input_refs_to_genome_refs(input_refs, self.kbase_api)
        print('collected refs:', genome_refs_all)
        genome_to_ref = validate_genome_refs(genome_refs_all)
        print('Genome To Ref:', genome_to_ref)

        messages = []
        warnings = []

        # Validate required parameters
        self._validate_params(params, ['input_refs', 'workspace_name'])

        with open(str(input_params.resolve()), 'w') as fh:
            _params = dict(params)
            _params['_ctx'] = ctx
            _params['_config'] = self.config
            _params['_genome_refs'] = genome_refs_all
            if DEBUG_MODE:
                print('to create a copy for debug:', _params)
            fh.write(json.dumps(_params))

        if DEBUG_MODE:
            print('data dir')
            print(os.listdir('/data'))
            if os.path.exists('/data') and os.path.exists('/data/reference_data'):
                print(os.listdir('/data/reference_data'))

                if os.path.exists('/data/reference_data/berdl_db'):
                    print('berdl:', os.listdir('/data/reference_data/berdl_db'))

            print('BERDL Token')
            print(get_berdl_token())

        #if not skip_annotation:
        #    test_annotation(self.kb_kofam, self.kb_bakta, self.kb_psortb, self.rast_client)

        input_genome_to_clade = {}
        clade_to_input_genomes = {}
        genome_refs = {}
        genome_refs_skipped = {}
        path_root = Path(self.shared_folder)
        if not skip_genome_pipeline:
            self.run_genome_pipeline(input_params.resolve())
            path_user_to_clade_json = path_root / 'pangenome' / 'user_to_clade.json'
            if path_user_to_clade_json.exists():
                print(f'found input genome pangenome assignments: {path_user_to_clade_json}')
                with open(path_user_to_clade_json, 'r') as fh:
                    input_genome_to_clade = json.load(fh)
            else:
                print(f'{path_user_to_clade_json} NOT FOUND')

            for genome_name, clade in input_genome_to_clade.items():
                if clade not in clade_to_input_genomes:
                    if len(clade_to_input_genomes) < param_clade_limit:
                        clade_to_input_genomes[clade] = set()
                    else:
                        warnings.append(f'WARNING CLADE LIMIT REACHED. Genome -> Clade skip. {genome_name} -> {clade}')
                        print(f'WARNING CLADE LIMIT REACHED. Genome -> Clade skip. {genome_name} -> {clade}')
                if clade in clade_to_input_genomes:
                    messages.append(f'Clade Match: {genome_name} -> {clade}')
                    clade_to_input_genomes[clade].add(genome_name)
                    genome_refs[genome_to_ref[genome_name]] = genome_name
                else:
                    genome_refs_skipped[genome_to_ref[genome_name]] = genome_name

            for genome_ref in genome_refs:
                info = self.util.get_object_info(genome_ref)
                path_genome_tsv = path_root / "genome" / f'user_{info[1]}_genome.tsv'
                print(f'create genome tsv: {path_genome_tsv} for {genome_ref}')
                self.util.run_user_genome_to_tsv(genome_ref, str(path_genome_tsv))
        else:
            print('skip genome pipeline')

        executor = TaskExecutor(max_workers=4)
        path_user_genome = path_root / "genome"
        path_user_genome.mkdir(parents=True, exist_ok=True)
        tasks_input_genome = []
        tasks_rast = []
        genome_allowed = {'user_' + g for g in genome_refs.values()}
        for filename_faa in path_user_genome.glob("*.faa"):
            genome_name = filename_faa.stem
            print(f'found protein faa {genome_name} -> {filename_faa}')
            if skip_annotation:
                print('skip_annotation')
            else:
                if genome_name in genome_allowed:
                    messages.append(f'Run annotation RAST: {filename_faa}')
                    th = executor.run_task(task_rast, path_user_genome / filename_faa, self.rast_client)
                    tasks_rast.append(th)
                    tasks_input_genome.append(th)
                    messages.append(f'Run annotation KOFAM: {filename_faa}')
                    tasks_input_genome.append(executor.run_task(task_kofam,
                                                                path_user_genome / filename_faa,
                                                                self.kb_kofam))
                    messages.append(f'Run annotation BAKTA: {filename_faa}')
                    tasks_input_genome.append(executor.run_task(task_bakta,
                                                                path_user_genome / filename_faa,
                                                                self.kb_bakta))
                else:
                    warnings.append(f"Skip annotation for: {genome_name}")
                    print(f'skip genome: {genome_name}')

        path_pangenome = path_root / "pangenome"
        path_pangenome.mkdir(parents=True, exist_ok=True)
        for folder_pangenome in os.listdir(str(path_pangenome)):
            if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                print(f'Found pangenome folder: {folder_pangenome}')
                # run pangenome pipeline for - folder_pangenome
                if not skip_pangenome:
                    if folder_pangenome == 'none':
                        warnings.append('Input genome contains genome(s) with no pangenome match.')
                    self.run_pangenome_pipeline(input_params.resolve(), folder_pangenome)
                else:
                    print('skip pangenome')

        tasks_pangeome = []
        path_pangenome = path_root / "pangenome"
        path_pangenome.mkdir(parents=True, exist_ok=True)
        for folder_pangenome in os.listdir(str(path_pangenome)):
            if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                print(f'Found pangenome folder: {folder_pangenome}')
                path_pangenome_members = path_pangenome / folder_pangenome / 'genome'
                if path_pangenome_members.exists() and not skip_annotation:
                    for _f in os.listdir(str(path_pangenome_members)):
                        if _f.endswith('.faa'):
                            th = executor.run_task(task_rast, path_pangenome_members / _f, self.rast_client)
                            tasks_rast.append(th)
                            tasks_pangeome.append(th)
                            tasks_pangeome.append(executor.run_task(self.run_annotation_pipeline,
                                                                    path_pangenome_members / _f))

        print('Task barrier input genome annotation RAST')
        for t in tasks_rast:
            print(f'await for {t.args} {t.status}')
            t.wait()

        
        if not skip_annotation:
            template_classifier = get_classifier()
            print(f'loaded classifier {template_classifier}')
            if template_classifier is not None:
                input_genome_class = {}
                for filename_faa in path_user_genome.iterdir():
                    genome_name = filename_faa.stem
                    if str(filename_faa).endswith('.faa') and genome_name in genome_allowed:
                        print(f'found faa {filename_faa}')
                        filename_faa_rast = path_user_genome / (filename_faa.stem + '_rast.tsv')
                        print(f'looking for {filename_faa_rast}')
                        if filename_faa_rast.exists():
                            print(f'found {filename_faa} with RAST: {filename_faa_rast}')
                            genome = read_rast_as_genome(filename_faa_rast, None)
                            res = template_classifier.classify(genome)
                            input_genome_class[genome_name] = res
                            messages.append(f'Genome Type Class: {filename_faa} -> {res}')
                            print(filename_faa, res)
                            psortb_org_param = '-n'
                            if res == 'P':
                                psortb_org_param = '-p'
                            elif res == 'A':
                                psortb_org_param = '-a'
                            messages.append(f'Run annotation PSORTB {psortb_org_param}: {filename_faa}')
                            tasks_input_genome.append(executor.run_task(task_psortb,
                                                                        path_user_genome / filename_faa,
                                                                        psortb_org_param,
                                                                        self.kb_psortb))
                with open(path_user_genome / 'genome_class.tsv', 'w') as fh:
                    h = ["genome_ref", "genome_name", "modelseed_class"]
                    fh.write("\t".join(h) + '\n')
                    for genome_ref, genome_name in genome_refs.items():
                        line = [genome_ref, genome_name, input_genome_class.get(genome_name, '?')]
                        fh.write("\t".join(line) + '\n')

        if not skip_modeling_pipeline:
            model_params = {
                "input_refs": genome_refs,
                "token": ctx['token'],
                "scratch": str(self.shared_folder),
                "kbversion": self.util.kb_version,
                "max_phenotypes": None,
                "module_path": "/kb/module"
            }
            model_params_file = path_root / 'model_pipeline_params.json'
            with open(str(model_params_file), 'w') as f:
                json.dump(model_params, f, indent=2)
            print(f"Wrote model pipeline params: {model_params_file}")
            self.run_model_pipeline(str(model_params_file))
        else:
            print('skip modeling pipeline')

        # Done with all tasks
        print('Task barrier input genome annotation')
        for t in tasks_input_genome:
            print(f'await for {t.args} {t.status}')
            t.wait()
        for t in tasks_input_genome:
            print(t.status)
            print(t.result)
            print(t.traceback)

        print('Task barrier')
        for t in tasks_pangeome:
            print(f'await for {t.args} {t.status}')
            t.wait()
        for t in tasks_pangeome:
            print(t.status)
            print(t.result)
            print(t.traceback)

        executor.shutdown()

        # safe to build table all task barrier reached
        """
        print(f'Export tsv tables [models, phenotypes]')
        self.util.build_model_tables(model_path=str(path_root / 'models'))
        self.util.build_phenotype_tables(
            output_dir=str(path_root / 'phenotypes'),
            phenosim_directory=str(path_root / 'phenotypes'),
            experiment_data_file='/kb/module/data/experimental_data.json',

            fitness_mapping_dir=str(path_root / 'genome'),
            model_data_dir=str(path_root / 'models'),

            fitness_genomes_dir='/data/reference_data/phenotype_data',
            reference_phenosim_dir='/data/reference_data/phenotype_data/phenosims'
        )
        """
        for folder_pangenome in os.listdir(str(path_pangenome)):
            if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                print(f'Build table for pangenome folder: {folder_pangenome}')
                # run table assembly pipeline for - folder_pangenome
                self.run_build_table(input_params.resolve(), folder_pangenome)
                path_db_file = path_root / 'pangenome' / folder_pangenome / 'db.sqlite'
                generate_ontology_tables(str(path_db_file),
                                         reference_data_path='/data/',
                                         source_tables=[
                                             'user_feature', 'pangenome_feature'
                                         ])

        print_path(path_root.resolve())



        # Save annotated genomes back to workspace from each clade's database
        genome_set_items = []
        if not skip_save_genome_annotation:
            for folder_pangenome in os.listdir(str(path_pangenome)):
                if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                    path_db_file = path_root / 'pangenome' / folder_pangenome / 'db.sqlite'
                    input_genome_names = clade_to_input_genomes[folder_pangenome]
                    input_genome_refs = [genome_to_ref[x] for x in input_genome_names]
                    if path_db_file.exists():
                        print(f'Saving annotated genomes {input_genome_refs} from {folder_pangenome}')
                        saved_items = self.util.save_annotated_genomes(
                            genome_refs=input_genome_refs,
                            suffix=suffix,
                            output_workspace=workspace_name,
                            database_filename=str(path_db_file),
                        )
                        genome_set_items += saved_items['items']

        print(f'saved items: {genome_set_items}')
        if len(genome_set_items) == 0:
            genome_set_items = [{"ref": ref, "label": ref} for ref in genome_refs]
        ref_genome_set = self.util.save_genome_set(f"annotated_genomes_{suffix}",
                                                   genome_set_items, workspace_name)

        # Safe to read and export data
        file_links = []
        # Zip up shared_folder contents and upload to Shock for downloadable report link
        """
        shared_folder_path = Path(self.shared_folder)
        if export_all_data:
            self.logger.info(f"Zipping shared folder contents: {shared_folder_path}")
            archive_shock_info = self.dfu.file_to_shock({
                'file_path': str(shared_folder_path),
                'pack': 'zip'
            })
            archive_shock_id = archive_shock_info['shock_id']
            self.logger.info(f"Uploaded shared folder archive to Shock: {archive_shock_id}")
            file_links.append({
                'shock_id': archive_shock_id,
                'name': 'pipeline_output.zip',
                'label': 'Pipeline Output',
                'description': 'Zipped archive of all pipeline output files'
            })
        """
        if export_genome_data:
            _file_path = path_root / 'genome'
            if _file_path.exists():
                if any(_file_path.iterdir()):
                    archive_shock_id = self.dfu.file_to_shock({
                        'file_path': str(_file_path),
                        'pack': 'zip'
                    })['shock_id']
                    file_links.append({
                        'shock_id': archive_shock_id,
                        'name': 'input_genomes.zip',
                        'label': 'Input Genomes Data',
                        'description': 'Input Genomes with annotation and model files'
                    })
        if export_folder_models:
            _file_path = path_root / 'models'
            if _file_path.exists():
                if any(_file_path.iterdir()):
                    archive_shock_id = self.dfu.file_to_shock({
                        'file_path': str(_file_path),
                        'pack': 'zip'
                    })['shock_id']
                    file_links.append({
                        'shock_id': archive_shock_id,
                        'name': 'models.zip',
                        'label': 'Models Folder',
                        'description': 'debug'
                    })
        if export_folder_phenotypes:
            _file_path = path_root / 'phenotypes'
            if _file_path.exists():
                if any(_file_path.iterdir()):
                    archive_shock_id = self.dfu.file_to_shock({
                        'file_path': str(_file_path),
                        'pack': 'zip'
                    })['shock_id']
                    file_links.append({
                        'shock_id': archive_shock_id,
                        'name': 'phenotypes.zip',
                        'label': 'Phenotypes Folder',
                        'description': 'debug'
                    })

        if export_pangenome_data:
            for folder_pangenome in path_pangenome.iterdir():
                if folder_pangenome.is_dir():
                    path_mmseqs_tmp = folder_pangenome / 'master_mmseqs2' / 'mmseqs2_tmp'
                    if os.path.exists(path_mmseqs_tmp):
                        shutil.rmtree(path_mmseqs_tmp)
                    name = folder_pangenome.name
                    print(f'found db for {folder_pangenome}')
                    archive_shock_id = self.dfu.file_to_shock({
                        'file_path': str(folder_pangenome),
                        'pack': 'zip'
                    })['shock_id']
                    file_links.append({
                        'shock_id': archive_shock_id,
                        'name': f'{name}.zip',
                        'label': f'{name} pangenome folder',
                        'description': f'{name} clade folder'
                    })

        if export_databases:
            path_export = path_root / 'export'
            path_export.mkdir(exist_ok=True)
            to_shock = {}
            for folder_pangenome in os.listdir(str(path_pangenome)):
                if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                    path_export_db = path_export / folder_pangenome
                    path_export_db.mkdir(exist_ok=True)
                    path_filename_db = (path_pangenome / folder_pangenome / 'db.sqlite').resolve()
                    if path_filename_db.exists():
                        shutil.copy2(path_filename_db, path_export_db / path_filename_db.name)
                        to_shock[folder_pangenome] = path_export_db
            for folder_pangenome, path_export_db in to_shock.items():
                print(f'file_to_shock: {path_export_db}')
                archive_shock_id = self.dfu.file_to_shock({
                    'file_path': str(path_export_db),
                    'pack': 'zip'
                })['shock_id']
                file_links.append({
                    'shock_id': archive_shock_id,
                    'name': f'{folder_pangenome}.zip',
                    'label': f'{folder_pangenome} database',
                    'description': f'{folder_pangenome} clade database'
                })

        pangenome_data_list = []
        for folder_pangenome in os.listdir(str(path_pangenome)):
            if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                path_db = (path_pangenome / folder_pangenome / 'db.sqlite').resolve()
                if path_db.exists():
                    shock_id, handle_id = upload_blob_file(str(path_db),
                                                           ctx['token'],
                                                           self.config['shock-url'],
                                                           self.hs)
                    # read members.tsv
                    print(f'upload_blob_file {path_db}: {shock_id} {handle_id}')
                    pangenome_data = {
                        'pangenome_id': folder_pangenome,
                        'pangenome_taxonomy': '',
                        'user_genomes': [],
                        'datalake_genomes': [],
                        'sqllite_tables_handle_ref': handle_id
                    }
                    pangenome_data_list.append(pangenome_data)

        # Create KBaseFBA.GenomeDataLakeTables

        output_object = {
            'name': output_object_name,
            'description': '',
            'genomeset_ref': ref_genome_set,
            'pangenome_data': pangenome_data_list,
        }
        saved_object_info = self.kbase_api.save_object(output_object_name,
                                                       params['workspace_name'],
                                                       'KBaseFBA.GenomeDataLakeTables',
                                                       output_object,
                                                       meta={}
        )

        print("saved object to workspace", saved_object_info)

        # Create report with results
        report_client = KBaseReport(self.callback_url)

        output_directory = os.path.join(self.shared_folder, str(uuid.uuid4()))
        shutil.copytree('/kb/module/data/html', output_directory)

        # Write app-config.json so the DataTables Viewer knows which object to display
        # Use the first input reference as the UPA for the viewer
        app_config = {
            "upa": str(saved_object_info)
        }
        print("HTML app_config", app_config)
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

        html_links = [{
            'shock_id': shock_id,
            'name': 'index.html',
            'label': 'BERDL Tables',
            'description': 'BERDL Table Viewer'
        }]

        report_params = {
            'message': '\n'.join(messages),
            'warnings': warnings,
            'workspace_name': params['workspace_name'],
            'objects_created': [],
            'html_links': html_links,
            'direct_html_link_index': 0,
            'file_links': file_links,
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
