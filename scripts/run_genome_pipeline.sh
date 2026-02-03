#!/usr/bin/env bash
set -euo pipefail

# Activate the berdl_genomes environment
source /opt/env/berdl_genomes/bin/activate

# Run the genome pipeline (force env)
/opt/env/berdl_genomes/bin/python /kb/module/berdl/berdl/pipeline.py "$@"
