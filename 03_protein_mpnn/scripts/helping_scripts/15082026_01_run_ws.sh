#!/usr/bin/env bash
# Stage 03, script 1 (workstation branch): resolves ProteinMPNN + conda for
# this machine, then hands off to the shared pipeline body in
# 15082026_01_run_one_job.sh (sourced, not executed -- so STAGE, MPNN_ROOT,
# and this shell's already-activated conda env all carry over automatically).
#
# Usage:
#   ./helping_scripts/15082026_01_run_ws.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="$(cd "$SCRIPT_DIR/../.." && pwd)"
export STAGE

MPNN_ROOT="${MPNN_ROOT:-/home/karina/software/ProteinMPNN}"
if [ ! -f "$MPNN_ROOT/protein_mpnn_run.py" ]; then
    echo "ERROR: protein_mpnn_run.py not found at: $MPNN_ROOT/protein_mpnn_run.py"
    echo "       Re-run as: MPNN_ROOT=/actual/path/to/ProteinMPNN bash $(basename "${BASH_SOURCE[0]}")"
    exit 1
fi

# Confirmed via `conda env list` on fabioadmin-ISRMOT-273: a dedicated
# "proteinmpnn" env exists there. Override if it's not the right one:
# MPNN_CONDA_ENV=your_env_name bash 15082026_01_run_ws.sh
MPNN_CONDA_ENV="${MPNN_CONDA_ENV:-proteinmpnn}"
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u   # conda's own activate.d hooks (e.g. boltz_cuda_libs.sh) reference
         # unset vars -- incompatible with -u, same reason the cluster
         # script already guards its own conda activate this way.
conda activate "$MPNN_CONDA_ENV"
set -u   # restore nounset now that conda's hooks have finished running

source "$SCRIPT_DIR/15082026_01_run_one_job.sh"
