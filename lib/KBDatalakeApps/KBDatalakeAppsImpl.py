# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
import shutil

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil

from cobrakbase import KBaseAPI

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
    GIT_COMMIT_HASH = "38ab5c13d23b7876d213dbc747d958806a775883"

    #BEGIN_CLASS_HEADER
    def _validate_params(self, params, required_keys):
        """Validate that required parameters are present."""
        for key in required_keys:
            if key not in params or params[key] is None:
                raise ValueError(f"Required parameter '{key}' is missing")
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

        # Initialize KBUtilLib utilities
        self.dfu = DataFileUtil(self.callback_url)
        self.kbase_api = KBaseAPI(os.environ['KB_AUTH_TOKEN'], config=config)
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

        # Validate required parameters
        self._validate_params(params, ['input_refs', 'workspace_name'])

        workspace_name = params['workspace_name']
        input_refs = params['input_refs']
        suffix = params.get('suffix', '')
        save_models = params.get('save_models', 0)

        print(params)

        # Process the input references
        results_text = f"Building genome datalake tables.\n"
        results_text += f"Input references: {len(input_refs)} object(s)\n"
        results_text += f"Suffix: {suffix if suffix else '(none)'}\n"
        results_text += f"Save models: {'Yes' if save_models else 'No'}\n\n"

        # TODO: Implement actual table building logic here
        # For now, just list the input refs
        results_text += "Input objects:\n"
        for ref in input_refs:
            results_text += f"  - {ref}\n"

        # Create KBaseFBA.GenomeDataLakeTables

        output_object = {
            'name': '',
            'description': '',
            'genomeset_ref': '',
            'pangenome_data': [{
                'pangenome_id': 'none',
                'pangenome_taxonomy': 'none',
                'user_genomes': [],
                'datalake_genomes': [],
                'table_handle_ref': ''
            }],
        }

        self.kbase_api.save_object('demo_output',
                                   params['workspace_name'],
                                   'KBaseFBA.GenomeDataLakeTables',
                                   output_object,
                                   meta=output_object # adding all data as metadata
        )

        # Create report with results
        report_client = KBaseReport(self.callback_url)

        output_directory = os.path.join(self.shared_folder, str(uuid.uuid4()))
        shutil.copytree('/kb/module/data/html', output_directory)

        shock_id = self.dfu.file_to_shock({
            'file_path': output_directory,
            'pack': 'zip'
        })['shock_id']

        print(os.listdir('/kb/module/data/html'))
        print(output_directory)
        print(os.listdir(output_directory))
        print(shock_id)
        

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
