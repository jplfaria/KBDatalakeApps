# -*- coding: utf-8 -*-
#BEGIN_HEADER
import json
import logging
import os
import uuid
from pathlib import Path
import shutil
import json
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
from modelseedpy import MSGenome
from cobrakbase import KBaseAPI
from installed_clients.baseclient import ServerError
from annotation.annotation import test_annotation

# Import KBUtilLib utilities for common functionality
#from kbutillib import KBWSUtils, KBCallbackUtils, SharedEnvUtils

#class DatalakeAppUtils(KBWSUtils, KBCallbackUtils, SharedEnvUtils):
#    """Custom utility class combining KBUtilLib modules for datalake operations."""
#    pass
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
    def get_berdl_token():
        return os.environ.get('KBASE_SECURE_CONFIG_PARAM_kbaselakehouseserviceaccount_token')

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

        print('polars thread pool', pl.thread_pool_size())
        self.rast_client = RAST_SDK(self.callback_url, service_ver='beta')
        #self.utils = DatalakeAppUtils(callback_url=self.callback_url)
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

        test_annotation(self.kb_kofam, self.kb_bakta, self.kb_psortb, self.rast_client)
        """
        def test_annotation_rast():
            # Printing test file for RAST annotation demonstration
            proteins = [
                ("Test3.CDS.1", "tRNA:Cm32/Um32 methyltransferase", "LFILTATGNMSLCGLKKECLIAASELVTCRE"),
                ("Test3.CDS.2", "Aspartokinase (EC 2.7.2.4);Homoserine dehydrogenase (EC 1.1.1.3)",
                 "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQEELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNKKANTSSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDEGMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLSWKLGV"),
            ]
            with open(self.shared_folder + "/test.faa", "w") as f:
                for seq_id, function, sequence in proteins:
                    f.write(f">{seq_id} {function}\n{sequence}\n")

            with open(self.shared_folder + "/test.faa", 'r') as fh:
                print('example faa:\n', fh.read())
            self.run_RAST_annotation(self.shared_folder + "/test.faa", self.shared_folder + "/rast.tsv",
                                     self.rast_client)

        try:
            test_annotation_rast()
        except ServerError as ex_server:
            logging.warning(f'error: {ex_server}')
        """

        #print('BERDL Token')
        #print(self.get_berdl_token())

        # Validate required parameters
        self._validate_params(params, ['input_refs', 'workspace_name'])

        workspace_name = params['workspace_name']
        input_refs = params['input_refs']
        suffix = params.get('suffix', ctx['token'])
        save_models = params.get('save_models', 0)

        self.run_genome_pipeline(input_params.resolve())

        path_pangenome = Path(self.shared_folder) / "pangenome"
        for folder_pangenome in os.listdir(str(path_pangenome)):
            if os.path.isdir(f'{path_pangenome}/{folder_pangenome}'):
                print(f'Found pangenome folder: {folder_pangenome}')
                # run pangenome pipeline for - folder_pangenome
                if not skip_pangenome:
                    self.run_pangenome_pipeline(input_params.resolve(), folder_pangenome)
                else:
                    print('skip pangenome')

        path_user_genome = Path(self.shared_folder) / "genome"
        t_start_time = time.perf_counter()
        for filename_faa in os.listdir(str(path_user_genome)):
            print('found', filename_faa)
            if skip_annotation:
                print('skip_annotation')
            else:
                if filename_faa.endswith('.faa'):
                    genome = MSGenome.from_fasta(str(path_user_genome / filename_faa))
                    proteins = {f.id:f.seq for f in genome.features if f.seq}
                    print('running annotation for', filename_faa, len(proteins))

                    try:
                        print(f"run kb_kofam annotation for {genome}")
                        self.logger.info(f"run annotation for {genome}")
                        start_time = time.perf_counter()
                        result = self.kb_kofam.annotate_proteins(proteins)
                        end_time = time.perf_counter()
                        print(f"Execution time: {end_time - start_time} seconds")
                        print(f'received results of type {type(result)} and size {len(result)}')

                        print(result)
                    except Exception as ex:
                        print(f'nope {ex}')

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

        t_end_time = time.perf_counter()
        print(f"Total Execution time annotation: {t_end_time - t_start_time} seconds")

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

        saved_object_info = self.kbase_api.save_object('fake_output',
                                   params['workspace_name'],
                                   'KBaseFBA.GenomeDataLakeTables',
                                   output_object,
                                   meta={
                                       'key1': 'value1',
                                       'note': 'all fields are fake values for testing'
                                   }
        )

        print(saved_object_info)

        # Create report with results
        report_client = KBaseReport(self.callback_url)

        output_directory = os.path.join(self.shared_folder, str(uuid.uuid4()))
        shutil.copytree('/kb/module/data/html', output_directory)

        # Write app-config.json so the DataTables Viewer knows which object to display
        # Use the first input reference as the UPA for the viewer
        app_config = {
            "upa": input_refs[0] if input_refs else None
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
