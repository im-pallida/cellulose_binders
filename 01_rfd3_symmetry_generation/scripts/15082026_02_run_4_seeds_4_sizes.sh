#!/usr/bin/env bash

set -euo pipefail

STAGE="$(cd "$(dirname"$BASH_SOURCE[0]")../" && pwd)"
MAX_CONCURRENT=3 # what is this value for?
ALL_SIZES=(small medium large xlarge)
ALL_ORDERS=(th_str ht_inv sh_str hs_inv)

EXPERIMENT="${1:-}" # means use the first value or we have nothing by default

if [ -z "$EXPERIMENT" ]; then # checking if the experiment is zero, if so
    "Usage: $0 <experiment>" # return this message
    exit 1 # means we exit bc of the general/miscellaneous error
fi
shift
MODE="auto" # what is this value for?
for order in "${ALL_ORDERS[@]}"; do
    if [ ! -d "$STAGE/json/$EXPERIMENT/$order" ]; then # if not statement, checking for the presence of jsons
        echo "ERROR: no such experiment/order folder: $STAGE/json/$EXPERIMENT/$order"
        echo "  (run the JSON-generation script first to create it)"
        exit 1
    fi
done
declare -A COUNT # declare an associative array
for order in "${ALL_ORDERS[@]}"; do
    for size in "${ALL_SIZES[@]}"; do
        read -r -p "How many '$order/$size' structures do you want to generate? (0 to skip): " val
        COUNT["${order}_${size}"]="$val"
        if ! [[ "${COUNT[${order}_${size}]}" =~ ^[0-9]+$ ]]; then
            echo "ERROR: count for '$order/$size' must be a non-negative integer, got: '${COUNT[${order}_${size}]}'"
            exit 1
        fi
    done
done
GLOBAL_COUNTER_FILE="$STAGE/job_runs/global_seq_counter.txt" # we count and record global sequences
mkdir -p "$(dirname "$GLOBAL_COUNTER_FILE")"
if [ -f "$GLOBAL_COUNTER_FILE" ]; then
    global_seq=$(cat "$GLOBAL_COUNTER_FILE")
else
    global_seq=0
fi
RUN_LIST="$STAGE/job_runs/${EXPERIMENT}_run.tsv"
: > "$RUN_LIST"
TOTAL=0
for order in "${ALL_ORDERS[@]}"; do
    for size in "${ALL_SIZES[@]}"; do
        count="${COUNT[${order}_${size}]}"
        if [ "$count" -eq 0 ]; then
            echo "[plan] $order/$size: skipped"
            continue
        fi
        start=$((global_seq + 1))
        for ((i = 0; i < count; i++)); do
            global_seq=$((global_seq + 1))
            printf "%s\t%s\t%d\n" "$order" "$size" "$global_seq" >> "$RUN_LIST"
            TOTAL=$((TOTAL + 1))
        done
        echo "[plan] $order/$size: $count structure(s), global seq $start-$global_seq"
    done
done
if [ "$TOTAL" -eq 0 ]; then
    echo "Nothing to do - all counts were 0."
    exit 0
fi
echo "$global_seq" > "$GLOBAL_COUNTER_FILE"
if [ "$MODE" = "auto" ]; then
    if command -v sbatch >/dev/null 2>&1; then
        MODE="--cluster"
    else
        MODE="--workstation"
    fi
fi
declare -A SIZE_CODE=( [small]="s" [medium]="m" [large]="l" [xlarge]="x" )
declare -A ORDER_CODE=( [th_str]="th" [ht_inv]="ht" [sh_str]="sh" [hs_inv]="hs" )
case "$MODE" in
    --workstation)
        echo "[dispatch] mode=workstation -> running $TOTAL job(s) locally, serial"
        exec "$STAGE/scripts/helping_scripts/15082026_02_run_ws.sh" "$EXPERIMENT" "$RUN_LIST"
        ;;
    --cluster)
        echo "[dispatch] mode=cluster -> canary job first, then $((TOTAL - 1)) more, max $MAX_CONCURRENT queued/running at once"
        # ---- CANARY: run job #1 with --wait, verify it actually succeeded ----
        read -r canary_order canary_size canary_seq < <(head -n 1 "$RUN_LIST")
        canary_ocode="${ORDER_CODE[$canary_order]}"
        canary_scode="${SIZE_CODE[$canary_size]}"
        canary_seq_padded=$(printf "%06d" "$canary_seq")
        canary_name="${canary_ocode}${canary_scode}${canary_seq_padded}"
        canary_raw_cif="$STAGE/outputs_raw/$EXPERIMENT/$canary_order/$canary_size/$canary_name.cif"
        echo "[canary] submitting job #1 ($canary_order $canary_size $canary_seq, name=$canary_name) with --wait ..."
        sbatch --wait "$STAGE/scripts/helping_scripts/15082026_02_run_cluster.sbatch" \
            "$EXPERIMENT" "$canary_order" "$canary_size" "$canary_seq"
        canary_exit=$?
        if [ "$canary_exit" -ne 0 ] || [ ! -f "$canary_raw_cif" ]; then
            echo
            echo "===== CANARY JOB FAILED - ABORTING, NOT SUBMITTING THE REMAINING $((TOTAL - 1)) JOB(S) ====="
            echo "  sbatch exit code: $canary_exit"
            echo "  expected output:  $canary_raw_cif"
            echo "  (this file does not exist - the job did not actually succeed, regardless of exit code)"
            echo "  Check the job's log under $STAGE/logs/ for the real error before resubmitting."
            echo "  Run 22072026_09_preflight_check.py $EXPERIMENT first to catch config-level issues (bad seed"
            echo "  paths, malformed json, etc.) before they cost you a whole batch again."
            exit 1
        fi
        echo "[canary] OK - job #1 succeeded, releasing the remaining $((TOTAL - 1)) job(s) at full concurrency"
        echo
        # ---- Remaining jobs: same as before, fire-and-poll at MAX_CONCURRENT ----
        tail -n +2 "$RUN_LIST" | while IFS=$'\t' read -r order size seq; do
            [ -n "$order" ] || continue
            while true; do
                in_queue=$(squeue -u "${USER:-$(whoami)}" -h | wc -l)
                if [ "$in_queue" -lt "$MAX_CONCURRENT" ]; then
                    break
                fi
                sleep 30
            done
            echo "[submit] $order $size $seq"
            sbatch "$STAGE/scripts/helping_scripts/15082026_02_run_cluster.sbatch" "$EXPERIMENT" "$order" "$size" "$seq"
            sleep 2
        done
        echo "All $TOTAL job(s) submitted (1 canary + $((TOTAL - 1)) batch)."
        ;;
    *)
        echo "ERROR: unknown mode '$MODE'"
        exit 1
        ;;
esac
