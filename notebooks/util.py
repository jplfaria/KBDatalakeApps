import sys
import os
import json
import shutil
from os import path

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
base_dir = os.path.dirname(os.path.dirname(script_dir))
folder_name = os.path.basename(script_dir)

print(base_dir+"/KBUtilLib/src")
sys.path = [base_dir+"/KBUtilLib/src",base_dir+"/cobrakbase",base_dir+"/ModelSEEDpy/"] + sys.path

# Add lib directory so KBDatalakeApps.KBDatalakeUtils can be imported
sys.path.insert(0, os.path.join(base_dir, "KBDatalakeApps", "lib"))

# Import utilities with error handling
from kbutillib import NotebookUtils, SharedEnvUtils, KBWSUtils

import pandas as pd


class NotebookUtil(NotebookUtils, KBWSUtils):
    """Test harness for the KBDatalakeApps pipeline.

    This utility reads token and config from the standard file locations
    (or from explicit arguments), and provides a helper to create
    KBDataLakeUtils instances with those credentials injected. This mirrors
    how the KBase SDK impl file initializes KBDataLakeUtils at runtime.
    """

    def __init__(self, **kwargs):
        super().__init__(
            notebook_folder=script_dir,
            name="KBDatalakeApps",
            user="chenry",
            retries=5,
            proxy_port=None,
            kb_version="appdev",
            **kwargs
        )

    def create_pipeline_utils(self, directory, parameters, kb_version, worker_count=4):
        """Create a KBDataLakeUtils instance with token/config injected from this notebook.

        This mirrors what the KBase SDK impl file does: it has access to the
        auth token and service config, and passes them into KBDataLakeUtils.
        The created instance will NOT read token/config from files itself.

        Args:
            directory: Working directory for pipeline output
            workspace_name: KBase workspace name for saving objects
            parameters: Pipeline parameters dict (must include 'input_refs')
            worker_count: Number of parallel workers (default: 4)

        Returns:
            Configured KBDataLakeUtils instance ready to run pipeline steps
        """
        from KBDatalakeApps.KBDatalakeUtils import KBDataLakeUtils

        # Pass all tokens (not just kbase) and the full config dict
        token = dict(self._token_hash) if self._token_hash else None
        config = dict(self._config_hash) if self._config_hash else None

        pipeline = KBDataLakeUtils(
            directory=directory,
            worker_count=worker_count,
            parameters=parameters,
            reference_path="/Users/chenry/Dropbox/Projects/KBDatalakeApps/notebooks/data/reference_data",
            # Inject token/config directly, skip all file reads
            token=token,
            config=config,
            token_file=None,
            kbase_token_file=None,
            config_file=False,
            kb_version=kb_version
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


# Initialize the NotebookUtil instance
util = NotebookUtil()
