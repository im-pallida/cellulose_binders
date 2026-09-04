#!/usr/bin/env bash
# Stage 04, script 2 (shared repetitive body): discovers pending AF3 jsons
# under $STAGE/json and, for each one not already done, runs AlphaFold3
# and records the result under $STAGE/outputs/<job_name>/. Sourced (not
# executed) by both 15082026_02_run_cluster.sbatch and
# 15082026_02_run_ws.sh -- so $STAGE and $AF3_EXEC_MODE ("cluster" or
# "workstation") must already be exported by the caller, along with
# whichever exec-mode-specific env vars that branch needs below
# (AF3_MODEL_PARAMETERS_DIR/AF3_PUB_DB_DIR/AF3_IMAGE_DIR/AF3_IMAGE for
# cluster, AF3_ROOT/MODEL_DIR for workstation).
#
# Idempotent: a job whose expected output already exists is skipped, same
# convention as every other *_run_one_job* script in this project -- so
# it's safe to rerun, or (on the cluster) to have several concurrent
# copies of this same loop running as separate pool workers.
#
# Concurrency safety (cluster only): before processing a json, a worker
# atomically claims it by `mkdir`-ing a marker under
# $STAGE/job_runs/.claims/<job_name> (mkdir fails if the dir already
# exists, so two workers racing on the same json can't both "win"). The
# claim is released again once that job is done (success or failure), so
# a failed job can still be retried on a later run. The workstation
# branch skips claiming entirely -- it's a single sequential process,
# there's nothing to race against.
# CAVEAT: if a worker is killed outright (node failure, OOM-kill) mid-job,
# its claim marker is left behind and that json won't be picked up again
# until you manually remove $STAGE/job_runs/.claims/<job_name>.
: "${STAGE:?STAGE not set -- this file must be sourced by run_cluster.sbatch or run_ws.sh, not run directly}"
: "${AF3_EXEC_MODE:?AF3_EXEC_MODE not set -- expected 'cluster' or 'workstation'}"

JSON_DIR="$STAGE/json"
OUTPUTS_DIR="$STAGE/outputs"
CLAIMS_DIR="$STAGE/job_runs/.claims"
mkdir -p "$OUTPUTS_DIR" "$CLAIMS_DIR"

echo "json      = $JSON_DIR"
echo "outputs   = $OUTPUTS_DIR"
echo "exec mode = $AF3_EXEC_MODE"

if [ ! -d "$JSON_DIR" ] || [ -z "$(ls -A "$JSON_DIR" 2>/dev/null)" ]; then
    echo "ERROR: $JSON_DIR is empty or missing -- run 15082026_01_af3_json_builder.py first"
    exit 1
fi

expected_result() {
    echo "$OUTPUTS_DIR/$1/${1}_model.cif"
}

already_done() {
    [ -f "$(expected_result "$1")" ]
}

claim_job() {
    local job_name="$1"
    if [ "$AF3_EXEC_MODE" != "cluster" ]; then
        return 0   # workstation: single process, nothing to race
    fi
    mkdir "$CLAIMS_DIR/$job_name" 2>/dev/null
}

release_claim() {
    local job_name="$1"
    [ "$AF3_EXEC_MODE" = "cluster" ] && rmdir "$CLAIMS_DIR/$job_name" 2>/dev/null
    return 0
}

run_af3_cluster() {
    local job_json="$1"
    local job_name="$2"
    local run_dir af3_input af3_output produced_dir dest_dir
    run_dir="$STAGE/af3_runs/$job_name"
    af3_input="$run_dir/af_input"
    af3_output="$run_dir/af_output_${SLURM_JOB_ID:-$$}"
    mkdir -p "$af3_input" "$af3_output"
    cp "$job_json" "$af3_input/fold_input.json"

    if ! apptainer exec \
        --nv \
        --bind "$af3_input:/root/af_input" \
        --bind "$af3_output:/root/af_output" \
        --bind "$AF3_MODEL_PARAMETERS_DIR:/root/models" \
        --bind "$AF3_PUB_DB_DIR:/root/public_databases" \
        "$AF3_IMAGE_DIR/$AF3_IMAGE" \
            python /app/alphafold/run_alphafold.py \
            --json_path=/root/af_input/fold_input.json \
            --model_dir=/root/models \
            --db_dir=/root/public_databases \
            --output_dir=/root/af_output; then
        return 1
    fi

    # AF3 names its output subfolder after the json's internal "name"
    # field (= job_name here), not the filename it was copied to.
    produced_dir="$af3_output/$job_name"
    if [ ! -d "$produced_dir" ]; then
        echo "  [ERROR] expected AF3 output not found at: $produced_dir"
        find "$af3_output" -maxdepth 2
        return 1
    fi
    dest_dir="$OUTPUTS_DIR/$job_name"
    rm -rf "$dest_dir"
    mv "$produced_dir" "$dest_dir"
    echo "  moved: $produced_dir -> $dest_dir"
}

run_af3_workstation() {
    local job_json="$1"
    # --norun_data_pipeline matches your original workstation launcher's
    # choice -- the cluster branch instead runs the full data pipeline
    # (--db_dir=...) per your explicit "accuracy over speed" choice there.
    python3 "$AF3_ROOT/run_alphafold.py" \
        --json_path="$job_json" \
        --model_dir="$MODEL_DIR" \
        --output_dir="$OUTPUTS_DIR" \
        --norun_data_pipeline
}

total_seen=0
skipped=0
succeeded=0
failed=()

for job_json in "$JSON_DIR"/*.json; do
    [ -e "$job_json" ] || continue
    job_name="$(basename "$job_json" .json)"
    total_seen=$((total_seen + 1))

    if already_done "$job_name"; then
        skipped=$((skipped + 1))
        continue
    fi
    if ! claim_job "$job_name"; then
        skipped=$((skipped + 1))   # another concurrent worker already claimed it
        continue
    fi
    # re-check after claiming -- another worker could have finished it
    # between our first check and winning the claim
    if already_done "$job_name"; then
        release_claim "$job_name"
        skipped=$((skipped + 1))
        continue
    fi

    echo
    echo "===== [$job_name] ====="
    ok=1
    if [ "$AF3_EXEC_MODE" = "cluster" ]; then
        run_af3_cluster "$job_json" "$job_name" || ok=0
    else
        run_af3_workstation "$job_json" || ok=0
    fi

    if [ "$ok" -eq 1 ] && already_done "$job_name"; then
        succeeded=$((succeeded + 1))
    else
        echo "  [WARNING] $job_name: did not produce $(expected_result "$job_name")"
        failed+=("$job_name")
    fi
    release_claim "$job_name"
done

echo
echo "Worker done. $total_seen job(s) seen, $skipped already-done/claimed-elsewhere (skipped), $succeeded newly succeeded."
if [ "${#failed[@]}" -gt 0 ]; then
    echo "[WARNING] ${#failed[@]} job(s) did not produce output: ${failed[*]}"
fi
echo "  outputs: $OUTPUTS_DIR"
