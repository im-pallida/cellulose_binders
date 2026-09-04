#!/usr/bin/env bash
#SBATCH --job-name=mpnn_launcher
#SBATCH --partition=gpu_l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
# NOTE: --output/--error are deliberately NOT set here (SLURM can't expand
# $STAGE in a #SBATCH comment at submit time) -- the dispatcher
# (15082026_01_run_proteinmpnn.sh) passes them on the sbatch command line
# instead, so logs land under this project's own logs/ dir wherever it
# lives on disk, rather than a path hardcoded to a different project.
#
# Stage 03, script 1 (cluster branch): resolves ProteinMPNN + conda for
# this machine, then hands off to the shared pipeline body in
# 15082026_01_run_one_job.sh (sourced, not executed -- so STAGE, MPNN_ROOT,
# and this shell's already-activated conda env all carry over automatically).
# This file only differs from the workstation branch in environment setup
# (conda activation, GPU node checks, resource limits above).
#
# NOTE: MPNN_ROOT / MINIFORGE_CONDA_SH / MPNN_CONDA_ENV_PATH below are
# inferred from your old project's cluster launcher (same cluster/
# allocation SCWF00146, same software install location under
# .../karina/software/, which is normally shared across your projects
# rather than per-project) -- double check these actually resolve on this
# project's cluster before trusting a real run. Override any of them via
# environment variables if they're wrong.
#
# SLURM copies THIS FILE into a per-job spool directory before running it
# (e.g. /var/spool/slurmd/job<id>/slurm_script), so ${BASH_SOURCE[0]} does
# NOT point at this script's real location the way it would for a normal
# bash invocation -- resolving STAGE/helping_scripts/ that way here
# silently breaks (this is exactly the "No such file or directory" you'd
# get trying to source 15082026_01_run_one_job.sh from the spool dir
# instead). The dispatcher (15082026_01_run_proteinmpnn.sh) passes both
# explicitly via `sbatch --export` instead, since ITS OWN path resolution
# is reliable (only the .sbatch file itself gets spool-copied by SLURM,
# not the dispatcher) -- so this file now REQUIRES submission through the
# dispatcher rather than a bare `sbatch 15082026_01_run_cluster.sbatch`.
#
# Usage:
#   ../15082026_01_run_proteinmpnn.sh --cluster
set -euo pipefail
HELPING_SCRIPTS_DIR="${HELPING_SCRIPTS_DIR:?HELPING_SCRIPTS_DIR not set -- submit via ../15082026_01_run_proteinmpnn.sh --cluster, not a bare sbatch on this file}"
STAGE="${STAGE:?STAGE not set -- submit via ../15082026_01_run_proteinmpnn.sh --cluster, not a bare sbatch on this file}"
export STAGE

echo "===== SBATCH JOB INFO ====="
echo "job id      = ${SLURM_JOB_ID:-N/A}"
echo "node        = $(hostname)"
echo "gpu check:"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || echo "  [WARNING] nvidia-smi not available on this node"
echo

MPNN_ROOT="${MPNN_ROOT:-/shared/scratch/SCWF00146/karina/software/ProteinMPNN}"
if [ ! -f "$MPNN_ROOT/protein_mpnn_run.py" ]; then
    echo "ERROR: protein_mpnn_run.py not found at: $MPNN_ROOT/protein_mpnn_run.py"
    echo "       This path is wrong for this node/session. Find the real path with:"
    echo "         find \"\$HOME\" -maxdepth 4 -iname protein_mpnn_run.py 2>/dev/null"
    echo "       Then re-run as: MPNN_ROOT=/actual/path sbatch 15082026_01_run_cluster.sbatch"
    exit 1
fi

echo
echo "===== ENVIRONMENT CHECK ====="
MINIFORGE_CONDA_SH="${MINIFORGE_CONDA_SH:-/shared/scratch/SCWF00146/karina/software/miniforge3/etc/profile.d/conda.sh}"
if [ ! -f "$MINIFORGE_CONDA_SH" ]; then
    echo "ERROR: expected miniforge3 conda.sh not found at: $MINIFORGE_CONDA_SH"
    echo "       Re-run as: MINIFORGE_CONDA_SH=/actual/path/conda.sh sbatch 15082026_01_run_cluster.sbatch"
    exit 1
fi
echo "[OK] sourcing miniforge3's own conda.sh directly: $MINIFORGE_CONDA_SH"
set +u   # conda's own activate.d hooks reference unset vars; incompatible with -u
source "$MINIFORGE_CONDA_SH"
MPNN_CONDA_ENV_PATH="${MPNN_CONDA_ENV_PATH:-/shared/scratch/SCWF00146/karina/software/miniforge3/envs/SE3nv}"
if [ ! -d "$MPNN_CONDA_ENV_PATH" ]; then
    echo "ERROR: expected conda environment not found at: $MPNN_CONDA_ENV_PATH"
    echo "       Re-run as: MPNN_CONDA_ENV_PATH=/actual/path sbatch 15082026_01_run_cluster.sbatch"
    set -u
    exit 1
fi
conda activate "$MPNN_CONDA_ENV_PATH"
set -u   # restore nounset now that conda's hooks have finished running
echo "[CHECK] confirming numpy is importable in this environment:"
python3 -c "import numpy; print('  numpy', numpy.__version__, 'at', numpy.__file__)"

source "$HELPING_SCRIPTS_DIR/15082026_01_run_one_job.sh"
