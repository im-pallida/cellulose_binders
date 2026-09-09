#!/usr/bin/env bash
#SBATCH --job-name=rfd3_stage01
#SBATCH --partition=gpu_h200
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=04:00:00
# NOTE: no --output/--error/--chdir here on purpose. SLURM does not expand shell
# variables in #SBATCH directives, so any path written here would be one user's
# absolute path. run_cluster.py passes them on the command line instead, derived
# from --stage, and command-line options override directives.
set -euo pipefail

# STAGE from this script's own location, so it does not matter who launches it
# or from which directory.
STAGE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# run_cluster.py submits exactly three arguments.
EXPERIMENT="$1"
JSON_REL="$2"
GLOBAL_SEQ="$3"

# run_one_job.py measures the seed with gemmi/numpy. Default to whichever python
# is active (the submitting user's conda env); override with RFD3_PIPELINE_PYTHON.
PYTHON="${RFD3_PIPELINE_PYTHON:-$(command -v python3 || true)}"
if [ -z "$PYTHON" ] || ! "$PYTHON" -c 'import gemmi, numpy' 2>/dev/null; then
    echo "ERROR: '${PYTHON:-<none>}' cannot import gemmi/numpy." >&2
    echo "       Activate the right conda env in this script, or set RFD3_PIPELINE_PYTHON." >&2
    exit 1
fi

exec "$PYTHON" "$STAGE/scripts/helping_scripts/run_one_job.py" \
     "$EXPERIMENT" "$JSON_REL" "$GLOBAL_SEQ" --stage "$STAGE"
