#!/usr/bin/env bash
# Stage 04, script 2 (cluster branch): submits AF3 jobs to SLURM, one
# `sbatch` per input json (15082026_02_run_one_job.sbatch), respecting
# your account's confirmed cap of 4 concurrent jobs/user (a 5th
# concurrent submission is rejected with QOSMaxSubmitJobPerUserLimit).
# Same throttled-submit logic as your working reference
# (28072026_02_submit_af3_jobs.sh) -- just split into functions here, and
# automated into one continuous run instead of two manual steps:
#   1. CANARY -- if nothing has ever completed successfully yet for this
#      stage, submit ONE job first (the first pending json) and wait for
#      it to actually produce output before touching anything else.
#   2. Once proven -- either just now, or already proven on a prior run --
#      queue ALL remaining pending jsons, throttled to stay at/under
#      MAX_CONCURRENT concurrent jobs: submits the next pending json
#      whenever squeue shows a free slot.
#
# Runs directly on the login node (bash, NOT itself an sbatch job) -- it
# does no compute itself, just sbatch + squeue + sleep in a loop, so it's
# safe to leave running in the foreground of a tmux/screen session for as
# long as it takes.
#
# Idempotent/resumable: a json already recorded in
# job_runs/af3_submitted.log, OR whose expected output already exists on
# disk (covers log drift or a manual resubmission -- same reasoning as
# the reference script), is skipped -- safe to Ctrl-C and rerun anytime.
#
# Usage:
#   ./helping_scripts/15082026_02_run_cluster.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="$(cd "$SCRIPT_DIR/../.." && pwd)"

JSON_DIR="$STAGE/json"          # written by 15082026_01_af3_json_builder.py
OUTPUTS_DIR="$STAGE/outputs"
LOGS_DIR="$STAGE/logs"
JOB_RUNS_DIR="$STAGE/job_runs"
SUBMITTED_LOG="$JOB_RUNS_DIR/af3_submitted.log"

# Found by pattern, not a hardcoded filename -- a renamed/date-drifted
# sbatch file caused this exact lookup to fail silently for hours in your
# reference project; matching by suffix survives that.
mapfile -t sbatch_candidates < <(find "$SCRIPT_DIR" -maxdepth 1 -name "*run_one_job.sbatch" | sort)
if [ "${#sbatch_candidates[@]}" -eq 0 ]; then
    echo "ERROR: no *run_one_job.sbatch found in $SCRIPT_DIR"
    exit 1
elif [ "${#sbatch_candidates[@]}" -gt 1 ]; then
    echo "ERROR: multiple *run_one_job.sbatch files found in $SCRIPT_DIR -- ambiguous:"
    printf '  %s\n' "${sbatch_candidates[@]}"
    exit 1
fi
SBATCH_SCRIPT="${sbatch_candidates[0]}"
echo "using sbatch template: $SBATCH_SCRIPT"

# Confirmed real limit on this cluster/allocation (SCWF00146) -- override
# only if you know it's changed.
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"
POLL_INTERVAL="${POLL_INTERVAL:-60}"

mkdir -p "$OUTPUTS_DIR" "$LOGS_DIR" "$JOB_RUNS_DIR"
touch "$SUBMITTED_LOG"

echo "json    = $JSON_DIR"
echo "outputs = $OUTPUTS_DIR"
echo "max concurrent = $MAX_CONCURRENT"

if [ ! -d "$JSON_DIR" ] || [ -z "$(ls -A "$JSON_DIR" 2>/dev/null)" ]; then
    echo "ERROR: $JSON_DIR is empty or missing -- run 15082026_01_af3_json_builder.py first"
    exit 1
fi

expected_result() {
    echo "$OUTPUTS_DIR/$1/${1}_model.cif"
}

any_completed_ever() {
    [ -n "$(find "$OUTPUTS_DIR" -mindepth 2 -maxdepth 2 -name "*_model.cif" -print -quit 2>/dev/null)" ]
}

already_submitted() {
    local job_json="$1"
    local job_name
    job_name="$(basename "$job_json" .json)"
    if grep -qxF "$job_json" "$SUBMITTED_LOG" 2>/dev/null; then
        return 0
    fi
    # Also treat a job as done if its real output already exists, even if
    # it's missing from the log (covers drift/duplicate-submission cases,
    # same reasoning as the reference script's already_submitted()).
    if [ -f "$(expected_result "$job_name")" ]; then
        echo "$job_json" >> "$SUBMITTED_LOG"
        return 0
    fi
    return 1
}

LAST_JOB_ID=""
submit_one() {
    local job_json="$1"
    local job_name sbatch_output
    job_name="$(basename "$job_json" .json)"
    if sbatch_output="$(sbatch \
        --job-name="af3_${job_name}" \
        --output="$LOGS_DIR/af3_${job_name}_%j.out" \
        --error="$LOGS_DIR/af3_${job_name}_%j.err" \
        "$SBATCH_SCRIPT" "$job_json" "$STAGE" 2>&1)"; then
        LAST_JOB_ID="$(echo "$sbatch_output" | grep -oE '[0-9]+$')"
        echo "  [submit] $job_name -> job $LAST_JOB_ID"
        echo "$job_json" >> "$SUBMITTED_LOG"
        return 0
    else
        echo "  [ERROR submitting] $job_name: $sbatch_output"
        return 1
    fi
}

current_job_count() {
    local output
    if ! output="$(squeue -u "$USER" -h 2>/dev/null)"; then
        echo "  [WARNING] squeue failed transiently -- assuming at limit, will retry" >&2
        echo "999999"   # force the caller to wait and re-check rather than crash
        return 0
    fi
    echo "$output" | wc -l
}

wait_for_success_or_failure() {
    # Waits until either any_completed_ever() is true, or job_id is no
    # longer in squeue -- re-checking any_completed_ever once more right
    # when the job disappears, since output landing on disk and squeue
    # dropping the job aren't strictly ordered from here.
    local job_id="$1"
    while true; do
        any_completed_ever && return 0
        if [ -z "$(squeue -j "$job_id" -h 2>/dev/null)" ]; then
            any_completed_ever && return 0
            return 1
        fi
        sleep "$POLL_INTERVAL"
    done
}

# Throttled submit: respects the cluster's per-user concurrent job limit.
# Submits one job at a time, waiting and re-polling squeue whenever
# at/above MAX_CONCURRENT, until every remaining json has been submitted.
submit_throttled() {
    local max_concurrent="$1"
    shift
    local jsons=("$@")
    local n_submitted=0
    local count
    echo "[plan] ${#jsons[@]} job(s) to submit, max $max_concurrent concurrent"
    for job_json in "${jsons[@]}"; do
        if already_submitted "$job_json"; then
            continue
        fi
        while true; do
            count="$(current_job_count)"
            if [ "$count" -lt "$max_concurrent" ]; then
                break
            fi
            sleep "$POLL_INTERVAL"
        done
        if submit_one "$job_json"; then
            n_submitted=$((n_submitted + 1))
        fi
        sleep 2
    done
    echo "[done] throttled submit finished: $n_submitted job(s) submitted this run"
}

# ---------------------------------------------------------------------------
# Main flow: canary (only if nothing has ever succeeded), then queue ALL
# remaining pending jsons, throttled.
# ---------------------------------------------------------------------------
mapfile -t all_jsons < <(find "$JSON_DIR" -maxdepth 1 -name "*.json" | sort)
echo "total input jsons found = ${#all_jsons[@]}"

pending=()
for job_json in "${all_jsons[@]}"; do
    already_submitted "$job_json" || pending+=("$job_json")
done
echo "pending job(s) to submit this run: ${#pending[@]}"
if [ "${#pending[@]}" -eq 0 ]; then
    echo "Nothing new to submit."
    exit 0
fi

if ! any_completed_ever; then
    canary_json="${pending[0]}"
    canary_name="$(basename "$canary_json" .json)"
    echo
    echo "===== CANARY: no successful AF3 result recorded yet -- submitting ONE job first: $canary_name ====="
    if ! submit_one "$canary_json"; then
        echo "ERROR: canary submission itself failed (sbatch error) -- aborting before submitting anything else."
        exit 1
    fi
    canary_job_id="$LAST_JOB_ID"
    echo "  waiting for job $canary_job_id to finish (polling every ${POLL_INTERVAL}s)..."
    if ! wait_for_success_or_failure "$canary_job_id"; then
        echo
        echo "ERROR: canary job $canary_name did not produce $(expected_result "$canary_name")."
        echo "       Check $LOGS_DIR/af3_${canary_name}_${canary_job_id}.err before rerunning this script."
        exit 1
    fi
    echo "  [OK] canary succeeded -- queueing the remaining $(( ${#pending[@]} - 1 )) job(s)."
    pending=("${pending[@]:1}")
fi

if [ "${#pending[@]}" -gt 0 ]; then
    echo
    submit_throttled "$MAX_CONCURRENT" "${pending[@]}"
fi
