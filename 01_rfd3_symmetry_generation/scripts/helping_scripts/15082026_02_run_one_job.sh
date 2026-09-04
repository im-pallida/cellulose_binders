#!/usr/bin/env bash
# Run exactly ONE structure. Given an experiment, an ORDER (th_str, ht_inv,
# sh_str, hs_inv), a size, and a GLOBAL sequence number (continuous across
# the whole project - not per-size, not per-order, not per-experiment),
# run RFD3, then:
#   - decompress + place the resulting structure at a flat path:
#       outputs_raw/<experiment>/<order>/<size>/<code><seq6>.cif
#     (plus its companion metadata file, same serial name:
#       outputs_raw/<experiment>/<order>/<size>/<code><seq6>.json)
#   - append both into a growing per-(experiment,order,size) archive:
#       outputs_clean/<experiment>/<order_code><size>.tar.gz
#     e.g. outputs_clean/<experiment>/htlarge.tar.gz
#
# RFD3 actually produces <example_id>_..._model_0.cif.gz (already gzipped)
# plus a sibling .json metadata file - confirmed against a real successful
# run, not guessed.
#
# ONLY .tar.gz is kept on disk (no uncompressed master .tar left behind).
# To append, an existing archive is decompressed to a throwaway temp .tar,
# extended with `tar -r`, re-gzipped, and the temp .tar is removed as part
# of that re-gzip step.
#
# Usage:
#   ./02_run_rfd3_one_job.sh <experiment> <order> <size> <global_seq>

set -euo pipefail

STAGE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXPERIMENT="$1"
ORDER="$2"
SIZE="$3"
GLOBAL_SEQ="$4"

declare -A SIZE_CODE=( [small]="s" [medium]="m" [large]="l" [xlarge]="x" )
declare -A ORDER_CODE=( [th_str]="th" [ht_inv]="ht" [sh_str]="sh" [hs_inv]="hs" )

SCODE="${SIZE_CODE[$SIZE]:-}"
if [ -z "$SCODE" ]; then
    echo "ERROR: unknown size '$SIZE' (expected small|medium|large|xlarge)"
    exit 1
fi

OCODE="${ORDER_CODE[$ORDER]:-}"
if [ -z "$OCODE" ]; then
    echo "ERROR: unknown order '$ORDER' (expected th_str|ht_inv|sh_str|hs_inv)"
    exit 1
fi

ENV_FILE="$STAGE/scripts/env/rfd3_t2exact.env"
OVERLAY_ROOT="$STAGE/overlay/rfd3_t3_overlay"

JSON_PATH="$STAGE/json/$EXPERIMENT/$ORDER/$SIZE.json"
SEQ_PADDED=$(printf "%06d" "$GLOBAL_SEQ")
NAME="${OCODE}${SCODE}${SEQ_PADDED}"

RAW_DIR="$STAGE/outputs_raw/$EXPERIMENT/$ORDER/$SIZE"
RAW_CIF="$RAW_DIR/$NAME.cif"
RAW_JSON="$RAW_DIR/$NAME.json"

CLEAN_DIR="$STAGE/outputs_clean/$EXPERIMENT"
CLEAN_ARCHIVE="$CLEAN_DIR/${OCODE}${SIZE}.tar.gz"
TMP_TAR="$CLEAN_DIR/.tmp_${OCODE}${SIZE}.tar"

TMP_WORK_DIR="$STAGE/outputs_raw/.tmp/$EXPERIMENT/$ORDER/$SIZE/$NAME"

if [ ! -f "$JSON_PATH" ]; then
    echo "ERROR: missing json config: $JSON_PATH"
    echo "  (did you run the JSON-generation script for this experiment?)"
    exit 1
fi

source "$ENV_FILE"

if [ -z "${RFD3_EXE:-}" ] || [ -z "${CKPT:-}" ]; then
    echo "ERROR: RFD3_EXE and/or CKPT not set - check $ENV_FILE defines both"
    exit 1
fi
if [ ! -x "$RFD3_EXE" ]; then
    echo "ERROR: RFD3_EXE does not exist or is not executable: $RFD3_EXE"
    exit 1
fi

export PYTHONPATH="$OVERLAY_ROOT:${PYTHONPATH:-}"

echo "===== RFD3 SIZES JOB ====="
date
hostname
echo "experiment=$EXPERIMENT  order=$ORDER  size=$SIZE  global_seq=$GLOBAL_SEQ  name=$NAME"
echo "json_path=$JSON_PATH"
echo "raw_cif=$RAW_CIF"
echo "clean_archive=$CLEAN_ARCHIVE"
echo "RFD3_T2EXACT_TRANSLATIONS=$RFD3_T2EXACT_TRANSLATIONS"
echo

if [ -f "$RAW_CIF" ]; then
    echo "[skip] $NAME -> $RAW_CIF already exists"
    exit 0
fi

echo "===== RUN RFD3 ====="
rm -rf "$TMP_WORK_DIR"
mkdir -p "$TMP_WORK_DIR"

"$RFD3_EXE" \
    inputs="$JSON_PATH" \
    out_dir="$TMP_WORK_DIR" \
    ckpt_path="$CKPT" \
    diffusion_batch_size=1 \
    n_batches=1 \
    skip_existing=False \
    inference_sampler.kind=symmetry \
    inference_sampler.num_timesteps=50 \
    +inference_sampler.sym_step_frac=0.9 \
    inference_sampler.use_classifier_free_guidance=False \
    inference_sampler.allow_realignment=False
echo
echo "===== LOCATE OUTPUT FILES ====="
mapfile -t cifgz_files < <(find "$TMP_WORK_DIR" -type f -iname "*.cif.gz")

if [ "${#cifgz_files[@]}" -eq 0 ]; then
    echo "ERROR: no .cif.gz file found in $TMP_WORK_DIR after RFD3 ran"
    find "$TMP_WORK_DIR" -type f
    exit 1
fi
if [ "${#cifgz_files[@]}" -gt 1 ]; then
    echo "ERROR: expected exactly 1 .cif.gz output, found ${#cifgz_files[@]}:"
    printf '  %s\n' "${cifgz_files[@]}"
    echo "Refusing to guess which one is correct."
    exit 1
fi

produced_cifgz="${cifgz_files[0]}"
echo "Found structure: $produced_cifgz"

# Don't assume the metadata json shares the exact basename with the
# .cif.gz - search the same directory for any .json file instead,
# mirroring the same robust find + count-check used for the .cif.gz above.
produced_dir="$(dirname "$produced_cifgz")"
mapfile -t json_files < <(find "$produced_dir" -maxdepth 1 -type f -iname "*.json")

json_found=false
if [ "${#json_files[@]}" -eq 0 ]; then
    echo "WARNING: no .json metadata file found in $produced_dir"
elif [ "${#json_files[@]}" -gt 1 ]; then
    echo "WARNING: expected exactly 1 .json metadata file in $produced_dir, found ${#json_files[@]}:"
    printf '  %s\n' "${json_files[@]}"
    echo "Refusing to guess which one is correct - metadata will NOT be archived."
else
    produced_json="${json_files[0]}"
    echo "Found metadata: $produced_json"
    json_found=true
fi

echo
echo "===== PLACE FLAT RAW FILES ====="
mkdir -p "$RAW_DIR"
gunzip -c "$produced_cifgz" > "$RAW_CIF"
echo "$RAW_CIF"

tar_members=("$NAME.cif")
if $json_found; then
    cp "$produced_json" "$RAW_JSON"
    echo "$RAW_JSON"
    tar_members+=("$NAME.json")
fi

echo
echo "===== APPEND TO GROWING ARCHIVE (tar.gz only) ====="
mkdir -p "$CLEAN_DIR"
rm -f "$TMP_TAR"

if [ -f "$CLEAN_ARCHIVE" ]; then
    echo "Extending existing archive: $CLEAN_ARCHIVE"
    gunzip -c "$CLEAN_ARCHIVE" > "$TMP_TAR"
    tar -rf "$TMP_TAR" -C "$RAW_DIR" "${tar_members[@]}"
else
    echo "Creating new archive: $CLEAN_ARCHIVE"
    tar -cf "$TMP_TAR" -C "$RAW_DIR" "${tar_members[@]}"
fi

gzip -f "$TMP_TAR"
mv -f "$TMP_TAR.gz" "$CLEAN_ARCHIVE"
ls -lh "$CLEAN_ARCHIVE"

echo
echo "===== CLEAN TEMP WORK DIR ====="
rm -rf "$TMP_WORK_DIR"
echo "[done] $NAME -> raw: $RAW_CIF, archive: $CLEAN_ARCHIVE"
date
