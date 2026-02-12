#!/usr/bin/env bash
set -euo pipefail

# Activate the berdl_genomes environment
source /opt/env/berdl_genomes/bin/activate
/opt/env/berdl_genomes/bin/python -c "import sys; print(sys.path)"

# Run the genome pipeline (force env)
/opt/env/berdl_genomes/bin/python /kb/module/berdl/berdl/bin/annotation.py "$@"
