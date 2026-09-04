#!/bin/bash --login
#SBATCH --job-name=af3_worker
#SBATCH --partition=gpu_l40s
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --mem=122G
#SBATCH --time=12:00:00
# NOTE: --output/--error are deliberately NOT set here (SLURM can't
# expand $STAGE in a #SBATCH comment at submit time) -- the dispatcher
# (15082026_02_run_alphafold.sh) passes them on the sbatch command line
# instead, same reason stage 03's cluster.sbatch does this too. Time
# limit is longer than a single-structure job (12h, not 4h) because this
# worker keeps looping and claiming more pending jsons until none are
# left, not just one.
#
# Stage 04, script 2 (cluster branch): ONE WORKER in a pool of up to 4
# concurrent copies (your account's confirmed cluster job cap). Each
# worker loops through pending AF3 jsons itself via the shared body in
# 15082026_02_run_one_job.sh (sourced, not executed), claiming one json
# at a time so concurrent workers never process the same one twice. This
# file only resolves the module/apptainer environment; the actual
# per-json AF3 invocation and claim logic lives in run_one_job.sh.
#
# STAGE / HELPING_SCRIPTS_DIR are required via `sbatch --export` from
# 15082026_02_run_alphafold.sh, NOT derived from ${BASH_SOURCE[0]} --
# SLURM copies this file into a per-job spool directory before running it
# (same spool-copy path-resolution issue already hit and fixed in stage
# 03's cluster.sbatch).
#
# Usage (normally only ever called by 15082026_02_run_alphafold.sh):
#   sbatch --export="ALL,STAGE=...,HELPING_SCRIPTS_DIR=..." 15082026_02_run_cluster.sbatch
set -euo pipefail
HELPING_SCRIPTS_DIR="${HELPING_SCRIPTS_DIR:?HELPING_SCRIPTS_DIR not set -- submit via ../15082026_02_run_alphafold.sh --cluster, not a bare sbatch on this file}"
STAGE="${STAGE:?STAGE not set -- submit via ../15082026_02_run_alphafold.sh --cluster, not a bare sbatch on this file}"
export STAGE

echo "===== SBATCH JOB INFO ====="
echo "job id = ${SLURM_JOB_ID:-N/A}"
echo "node   = $(hostname)"
echo "gpu check:"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || echo "  [WARNING] nvidia-smi not available on this node"
echo

module purge
module load alphafold/3.0.1
module list

# Inferred shared software location (same cluster/allocation SCWF00146,
# same "normally shared across your projects" software layout your
# ProteinMPNN cluster launcher already used) -- double check before
# trusting a real run; override via env var if wrong.
AF3_WORKDIR="${AF3_WORKDIR:-/scratch/SCWF00146/karina/software/alphafold}"
export AF3_MODEL_PARAMETERS_DIR="${AF3_MODEL_PARAMETERS_DIR:-$AF3_WORKDIR/model_parameters}"
# AF3_PUB_DB_DIR / AF3_IMAGE_DIR / AF3_IMAGE deliberately have NO default
# here -- same as your confirmed-working reference script, they're
# expected to already be set by `module load alphafold/3.0.1` (or your
# login profile). Failing fast beats apptainer silently binding an empty
# path if one of them isn't actually set.
: "${AF3_PUB_DB_DIR:?AF3_PUB_DB_DIR not set -- expected from 'module load alphafold/3.0.1'}"
: "${AF3_IMAGE_DIR:?AF3_IMAGE_DIR not set -- expected from 'module load alphafold/3.0.1'}"
: "${AF3_IMAGE:?AF3_IMAGE not set -- expected from 'module load alphafold/3.0.1'}"
export AF3_PUB_DB_DIR AF3_IMAGE_DIR AF3_IMAGE

export AF3_EXEC_MODE=cluster
source "$HELPING_SCRIPTS_DIR/15082026_02_run_one_job.sh"
