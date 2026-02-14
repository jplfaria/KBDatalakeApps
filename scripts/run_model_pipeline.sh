#!/usr/bin/env bash
set -euo pipefail

# Run model pipeline using SDK module Python (NOT berdl_genomes env)
python /kb/module/berdl/berdl/bin/model_pipeline.py "$@"
