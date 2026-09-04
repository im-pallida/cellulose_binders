#!/usr/bin/env bash
# Stage 04, script 2 (workstation branch): runs AF3 jobs one at a time,
# directly -- there's no SLURM queue/cap on the workstation, so no
# throttling or claim logic is needed here (unlike
# 15082026_02_run_cluster.sh's canary + max-4-concurrent submitter).
#
# Idempotent: a job whose expected output already exists is skipped, so
# this is safe to rerun after a partial/interrupted run.
#
# Usage:
#   ./helping_scripts/15082026_02_run_ws.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="$(cd "$SCRIPT_DIR/../.." && pwd)"

JSON_DIR="$STAGE/json"       # written by 15082026_01_af3_json_builder.py
OUTPUTS_DIR="$STAGE/outputs"
mkdir -p "$OUTPUTS_DIR"

echo "json    = $JSON_DIR"
echo "outputs = $OUTPUTS_DIR"

if [ ! -d "$JSON_DIR" ] || [ -z "$(ls -A "$JSON_DIR" 2>/dev/null)" ]; then
    echo "ERROR: $JSON_DIR is empty or missing -- run 15082026_01_af3_json_builder.py first"
    exit 1
fi

AF3_ROOT="${AF3_ROOT:-/home/karina/software/alphafold3}"
if [ ! -f "$AF3_ROOT/run_alphafold.py" ]; then
    echo "ERROR: run_alphafold.py not found at: $AF3_ROOT/run_alphafold.py"
    echo "       Re-run as: AF3_ROOT=/actual/path/to/alphafold3 bash $(basename "${BASH_SOURCE[0]}")"
    exit 1
fi
MODEL_DIR="${MODEL_DIR:-/home/fabioadmin/software/alphafold3_model_params}"
if [ ! -d "$MODEL_DIR" ]; then
    echo "ERROR: model params dir not found at: $MODEL_DIR"
    echo "       Re-run as: MODEL_DIR=/actual/path bash $(basename "${BASH_SOURCE[0]}")"
    exit 1
fi

AF3_CONDA_ENV="${AF3_CONDA_ENV:-alphafold3}"
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u   # conda's own activate.d hooks reference unset vars; incompatible with -u
conda activate "$AF3_CONDA_ENV"
set -u

expected_result() {
    echo "$OUTPUTS_DIR/$1/${1}_model.cif"
}

total=0
skipped=0
succeeded=0
failed=()

for job_json in "$JSON_DIR"/*.json; do
    [ -e "$job_json" ] || continue
    job_name="$(basename "$job_json" .json)"
    total=$((total + 1))
    result_path="$(expected_result "$job_name")"
    if [ -f "$result_path" ]; then
        skipped=$((skipped + 1))
        continue
    fi

    echo
    echo "===== [$job_name] ====="
    # --norun_data_pipeline matches your original workstation launcher's
    # choice -- the cluster branch instead runs the full data pipeline
    # (--db_dir=...) per your explicit "accuracy over speed" choice there.
    if python3 "$AF3_ROOT/run_alphafold.py" \
        --json_path="$job_json" \
        --model_dir="$MODEL_DIR" \
        --output_dir="$OUTPUTS_DIR" \
        --norun_data_pipeline; then
        if [ -f "$result_path" ]; then
            succeeded=$((succeeded + 1))
        else
            echo "  [WARNING] $job_name: run_alphafold.py exited 0 but $result_path is still missing"
            failed+=("$job_name")
        fi
    else
        echo "  [WARNING] $job_name: run_alphafold.py failed -- continuing to the next job"
        failed+=("$job_name")
    fi
done

echo
echo "Run complete. $total job(s) seen, $skipped already done (skipped), $succeeded newly succeeded."
if [ "${#failed[@]}" -gt 0 ]; then
    echo "[WARNING] ${#failed[@]} job(s) did not produce output: ${failed[*]}"
fi
echo "  outputs: $OUTPUTS_DIR"
