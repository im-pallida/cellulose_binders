#!/bin/bash --login
#SBATCH --job-name=af3_job
#SBATCH --partition=gpu_l40s
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --mem=122G
#SBATCH --time=04:00:00
# NOTE: --job-name/--output/--error are deliberately NOT relied on from
# these #SBATCH comments -- 15082026_02_run_cluster.sh passes them
# explicitly on the sbatch command line instead (SLURM can't expand
# $STAGE/$job_name inside a #SBATCH comment at submit time).
#
# The repetitive body: runs ONE AlphaFold3 job (one input json) and moves
# its result into <STAGE>/outputs/<job_name>/. This file itself doesn't
# loop over anything -- 15082026_02_run_cluster.sh submits it once per
# pending json, throttled to stay under your account's confirmed 4
# concurrent jobs/user. Adapted from your confirmed-working
# 28072026_02_af3_run_one_job.sbatch -- same module, same apptainer
# invocation, same --db_dir (full data pipeline, real MSA search, per
# your explicit accuracy-over-speed choice).
#
# JOB_JSON and STAGE are passed explicitly as positional args by
# 15082026_02_run_cluster.sh -- NOT derived from ${BASH_SOURCE[0]} here,
# because SLURM copies this file into a per-job spool directory before
# running it (the spool-copy path-resolution issue already hit and fixed
# in stage 03's cluster.sbatch -- this script's own location can't be
# trusted once SLURM has copied it).
#
# Usage (normally only ever called by 15082026_02_run_cluster.sh):
#   sbatch 15082026_02_run_one_job.sbatch /path/to/some_job.json /path/to/04_alphafold
set -eu
if [ $# -ne 2 ]; then
    echo "ERROR: usage: sbatch $0 /path/to/job.json /path/to/04_alphafold"
    exit 1
fi
JOB_JSON="$1"
STAGE="$2"
JOB_NAME="$(basename "$JOB_JSON" .json)"

FINAL_OUTPUTS_DIR="$STAGE/outputs"
EXPECTED_RESULT="$FINAL_OUTPUTS_DIR/$JOB_NAME/${JOB_NAME}_model.cif"

echo "job_name = $JOB_NAME"
echo "input    = $JOB_JSON"
echo "stage    = $STAGE"

# Same idempotency pattern as 16072026_02_run_rfd3_one_job.sh / your
# confirmed-working reference: check for the actual output before doing
# any real work, not just an external submission log -- protects against
# duplicate submissions from ANY cause (log drift, manual resubmission,
# restarted scripts) burning GPU time redoing a job that already
# succeeded.
if [ -f "$EXPECTED_RESULT" ]; then
    echo "[skip] $JOB_NAME -> $EXPECTED_RESULT already exists"
    exit 0
fi

module purge
module load alphafold/3.0.1
module list

# Inferred shared software location (same cluster/allocation SCWF00146,
# same "normally shared across your projects" software layout your
# ProteinMPNN cluster launcher already used) -- double check before
# trusting a real run; override via env var if wrong.
AF3_WORKDIR="${AF3_WORKDIR:-/scratch/SCWF00146/karina/software/alphafold}"
AF3_MODEL_PARAMETERS_DIR="${AF3_MODEL_PARAMETERS_DIR:-$AF3_WORKDIR/model_parameters}"
# AF3_PUB_DB_DIR / AF3_IMAGE_DIR / AF3_IMAGE deliberately have NO default
# here -- same as your confirmed-working reference script, they're
# expected to already be set by `module load alphafold/3.0.1` (or your
# login profile). Failing fast beats apptainer silently binding an empty
# path if one of them isn't actually set.
: "${AF3_PUB_DB_DIR:?AF3_PUB_DB_DIR not set -- expected from 'module load alphafold/3.0.1'}"
: "${AF3_IMAGE_DIR:?AF3_IMAGE_DIR not set -- expected from 'module load alphafold/3.0.1'}"
: "${AF3_IMAGE:?AF3_IMAGE not set -- expected from 'module load alphafold/3.0.1'}"

RUN_DIR="$STAGE/af3_runs/$JOB_NAME"
AF3_INPUT="$RUN_DIR/af_input"
AF3_OUTPUT="$RUN_DIR/af_output_${SLURM_JOB_ID}"
mkdir -p "$AF3_INPUT" "$AF3_OUTPUT" "$FINAL_OUTPUTS_DIR" "$STAGE/logs"
cp "$JOB_JSON" "$AF3_INPUT/fold_input.json"
echo "run_dir  = $RUN_DIR"

apptainer exec \
    --nv \
    --bind "$AF3_INPUT:/root/af_input" \
    --bind "$AF3_OUTPUT:/root/af_output" \
    --bind "$AF3_MODEL_PARAMETERS_DIR:/root/models" \
    --bind "$AF3_PUB_DB_DIR:/root/public_databases" \
    "$AF3_IMAGE_DIR/$AF3_IMAGE" \
        python /app/alphafold/run_alphafold.py \
        --json_path=/root/af_input/fold_input.json \
        --model_dir=/root/models \
        --db_dir=/root/public_databases \
        --output_dir=/root/af_output

echo
echo "===== CONSOLIDATE OUTPUT ====="
# AF3 names its output subfolder after the json's internal "name" field
# (= JOB_NAME here), not the filename it was copied to (fold_input.json).
PRODUCED_DIR="$AF3_OUTPUT/$JOB_NAME"
if [ ! -d "$PRODUCED_DIR" ]; then
    echo "ERROR: expected AF3 output not found at: $PRODUCED_DIR"
    echo "Listing what AF3 actually produced instead:"
    find "$AF3_OUTPUT" -maxdepth 2
    exit 1
fi
DEST_DIR="$FINAL_OUTPUTS_DIR/$JOB_NAME"
rm -rf "$DEST_DIR"
mv "$PRODUCED_DIR" "$DEST_DIR"
echo "moved: $PRODUCED_DIR -> $DEST_DIR"

echo
echo "[DONE] $JOB_NAME"
