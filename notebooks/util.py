import sys
import os
import json
import shutil
from os import path
from pathlib import Path

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
base_dir = "/home/chenry/Dropbox/Projects"
kbutillib_dir = os.path.dirname(script_dir)

sys.path = [base_dir+"/KBUtilLib/src",base_dir+"/cobrakbase",base_dir+"/ModelSEEDpy/","/home/chenry/MyEnvs/modelseed_cplex"] + sys.path

# Add lib directory so KBDatalakeApps.KBDatalakeUtils can be imported
sys.path.insert(0, os.path.join(kbutillib_dir,"lib"))

# Import utilities with error handling
from kbutillib import NotebookUtils, SharedEnvUtils,KBWSUtils

import pandas as pd

# Import CatalogClient from lib directory (already in sys.path from line above)
from kbase_catalog_client import CatalogClient, CatalogError

# Import the same functions from KBDatalakeUtils that the Impl file does:
#   from KBDatalakeApps.KBDatalakeUtils import KBDataLakeUtils, run_phenotype_simulation, run_model_reconstruction
from KBDatalakeApps.KBDatalakeUtils import KBDataLakeUtils, run_phenotype_simulation, run_model_reconstruction

# RAST annotation support (modelseedpy RastClient for notebook context,
# replaces the KBase RAST_SDK callback client used in the deployed app)
from modelseedpy.core.msgenome import MSGenome, MSFeature
from modelseedpy.core.rast_client import RastClient


class NotebookUtil(NotebookUtils,SharedEnvUtils):
    """Proxy for KBDatalakeAppsImpl for notebook-based testing.

    Replicates the Impl file's initialization of KBDataLakeUtils and imports
    the same functions. Acts as a test harness that mirrors the deployed
    KBase app as closely as possible.

    In the deployed app, the Impl constructor does:
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.config = config
        self.dfu = DataFileUtil(self.callback_url)
        self.kbase_api = KBaseAPI(os.environ['KB_AUTH_TOKEN'], config=config)
        self.rast_client = RAST_SDK(self.callback_url, service_ver='beta')
        self.util = KBDataLakeUtils(
            kbendpoint=config["kbase-endpoint"],
            reference_path="/data/",
            module_path="/kb/module"
        )
        self.util.set_token(get_berdl_token(), namespace="berdl")

    This notebook util replicates that pattern using local paths and
    file-based tokens instead of the SDK callback and deploy.cfg.
    """

    def __init__(self, **kwargs):
        super().__init__(
            notebook_folder=script_dir,
            name="KBDatalakeApps",
            **kwargs
        )

    def create_kbdl_utils(self, reference_path,module_path,kb_version):
        """Create a KBDataLakeUtils instance mirroring the Impl constructor.

        Impl does:
            self.util = KBDataLakeUtils(
                kbendpoint=config["kbase-endpoint"],
                reference_path="/data/",
                module_path="/kb/module"
            )
            self.util.set_token(get_berdl_token(), namespace="berdl")

        This method replicates that using local paths and notebook tokens.

        Args:
            reference_path: Override reference data path (default: notebooks/data/reference_data)
            module_path: Override module path (default: KBDatalakeApps root)

        Returns:
            Configured KBDataLakeUtils instance
        """
        utils = KBDataLakeUtils(
            reference_path=reference_path,
            module_path=module_path,
            kb_version=kb_version,
            token=self.get_token("kbase")
        )

        # Mirror: self.util.set_token(get_berdl_token(), namespace="berdl")
        berdl_token = self.get_token('berdl')
        if berdl_token:
            utils.set_token(berdl_token, namespace="berdl")

        return utils

    def create_pipeline_utils(self, directory, parameters, kb_version=None, worker_count=4):
        """Create a KBDataLakeUtils instance with directory and parameters injected.

        Used by test_pipeline_steps.ipynb for step-by-step pipeline testing.

        Args:
            directory: Working directory for pipeline output
            parameters: Pipeline parameters dict (must include 'input_refs')
            kb_version: KBase environment version (default: self.kb_version)
            worker_count: Number of parallel workers (default: 4)

        Returns:
            Configured KBDataLakeUtils instance
        """
        if kb_version is None:
            kb_version = self.kb_version

        token = dict(self._token_hash) if self._token_hash else None
        config = dict(self._config_hash) if self._config_hash else None

        pipeline = KBDataLakeUtils(
            kb_version=kb_version,
            reference_path=os.path.join(
                base_dir, "KBDatalakeApps", "notebooks", "data", "reference_data"
            ),
            module_path=os.path.join(base_dir, "KBDatalakeApps"),
            token=token,
            config=config,
            token_file=None,
            kbase_token_file=None,
            config_file=False,
        )

        return pipeline

    def clear_pipeline_dir(self, pipeline_dir):
        """Remove all files and subdirectories inside pipeline_dir.

        The directory itself is preserved so subsequent pipeline steps
        can write into it without needing to recreate it.

        Args:
            pipeline_dir: Path to the pipeline working directory to clear
        """
        if not os.path.isdir(pipeline_dir):
            self.log_warning(f"Pipeline directory does not exist: {pipeline_dir}")
            return
        for entry in os.listdir(pipeline_dir):
            entry_path = os.path.join(pipeline_dir, entry)
            if os.path.isdir(entry_path):
                shutil.rmtree(entry_path)
            else:
                os.remove(entry_path)
        self.log_info(f"Cleared pipeline directory: {pipeline_dir}")

    def get_catalog_client(self, environment=None):
        """Create a CatalogClient instance with the KBase token from this notebook.

        Args:
            environment: KBase environment (appdev, prod, ci, next).
                         Defaults to the notebook's kb_version setting.

        Returns:
            CatalogClient instance configured with KBase auth token
        """
        env = environment or getattr(self, 'kb_version', 'appdev')
        token = self.get_kbase_token()
        return CatalogClient(token=token, environment=env)

    def annotate_faa_with_rast(self, faa_path, output_tsv_path=None):
        """Annotate a FASTA protein file with RAST and write a 2-column TSV.

        Mirrors the Impl's run_RAST_annotation static method but uses
        modelseedpy's RastClient instead of the KBase RAST_SDK callback client.

        Args:
            faa_path: Path to input .faa file
            output_tsv_path: Path for output .tsv (default: same name with .tsv extension)

        Returns:
            Path to the output TSV file
        """
        faa_path = Path(faa_path)
        if output_tsv_path is None:
            output_tsv_path = faa_path.with_suffix('.tsv')
        else:
            output_tsv_path = Path(output_tsv_path)

        # Parse FASTA (same approach as Impl's run_RAST_annotation)
        sequences = {}
        current_id = None
        current_seq = []
        with open(faa_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        sequences[current_id] = "".join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                elif line:
                    current_seq.append(line)
            if current_id:
                sequences[current_id] = "".join(current_seq)

        # Create MSGenome and annotate with RAST
        genome = MSGenome()
        features = [MSFeature(fid, seq) for fid, seq in sequences.items()]
        genome.add_features(features)

        rast = RastClient()
        rast.annotate_genome(genome)

        # Write 2-column TSV: id, functions
        records = []
        for feature in genome.features:
            funcs = feature.ontology_terms.get('RAST', [])
            records.append({
                'id': feature.id,
                'functions': ';'.join(str(f) for f in funcs) if funcs else ''
            })

        df = pd.DataFrame(records)
        df.to_csv(output_tsv_path, sep='\t', index=False)
        return str(output_tsv_path)


# Initialize the NotebookUtil instance
util = NotebookUtil()
