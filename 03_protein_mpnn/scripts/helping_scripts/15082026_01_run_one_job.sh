#!/usr/bin/env bash
# Stage 03, script 1 (shared pipeline body): the actual ProteinMPNN run --
# identical between the workstation and cluster branches. Sourced (not
# executed) by 15082026_01_run_ws.sh / 15082026_01_run_cluster.sbatch
# AFTER they've each done their own machine-specific setup (STAGE export,
# MPNN_ROOT resolution + validation, conda activation). Being sourced
# means it runs in the caller's own shell: STAGE/MPNN_ROOT and the
# already-active conda env are inherited automatically, and the caller's
# `set -euo pipefail` stays in effect here too -- no need to redeclare it.
#
# Expects, already set by the caller before sourcing this file:
#   STAGE      -- 03_protein_mpnn/ root (exported)
#   MPNN_ROOT  -- validated ProteinMPNN install (protein_mpnn_run.py present)
# and an active conda env with ProteinMPNN's own dependencies.
#
# tied/symmetric design (chain A forced identical to chain B), pre-flight
# symmetry check on one structure before the full batch, cellulose chain
# stripped before parsing, output nested outputs/<experiment>/<group>.tar.gz
# (experiment looked up per protein_id from scripts/map/protein_experiment_map.csv,
# written by step 00 -- inputs_prepared/ itself is flat, so that file is
# the only place experiment provenance survives), cumulative archiving
# (decompress-append-recompress + atomic replace), raw scratch output
# cleaned up only after a confirmed-successful archive merge (step 5).

NUM_SEQ=16
INPUTS_PREPARED_DIR="$STAGE/inputs_prepared"
JSONLS_DIR="$STAGE/scripts/jsonls"
FIXED_POSITIONS_JSONL="$JSONLS_DIR/fixed_positions.jsonl"
BIAS_AA_JSONL="$JSONLS_DIR/bias_AA.jsonl"
TABLES_DIR="$STAGE/tables"
MAP_DIR="$STAGE/scripts/map"
PROTEIN_EXPERIMENT_MAP="$MAP_DIR/protein_experiment_map.csv"
export PROTEIN_EXPERIMENT_MAP
OUTPUTS_DIR="$STAGE/outputs"
RAW_SEQS_DIR="$STAGE/mpnn_raw"  # protein_mpnn_run.py's raw combined-fasta output, before splitting -- top level, not under tables/ (it's an output artifact, not a table)
MPNN_STAGING_DIR="$STAGE/.mpnn_staging"
export MPNN_STAGING_DIR
PENDING_COUNT_FILE="$STAGE/.mpnn_pending_count"
export PENDING_COUNT_FILE

echo "inputs_prepared  = $INPUTS_PREPARED_DIR"
echo "fixed_positions  = $FIXED_POSITIONS_JSONL"
echo "bias_AA          = $BIAS_AA_JSONL"
echo "experiment_map   = $PROTEIN_EXPERIMENT_MAP"
echo "outputs          = $OUTPUTS_DIR"
echo "mpnn_raw         = $RAW_SEQS_DIR (scratch -- cleared after each successful run)"

if [ ! -d "$INPUTS_PREPARED_DIR" ] || [ -z "$(ls -A "$INPUTS_PREPARED_DIR" 2>/dev/null)" ]; then
    echo "ERROR: $INPUTS_PREPARED_DIR is missing or empty -- run 15082026_00_pdb_extraction_jsonls_preparation.py first"
    exit 1
fi
if [ ! -f "$FIXED_POSITIONS_JSONL" ]; then
    echo "ERROR: $FIXED_POSITIONS_JSONL not found -- run 15082026_00_pdb_extraction_jsonls_preparation.py first"
    exit 1
fi
if [ ! -f "$BIAS_AA_JSONL" ]; then
    echo "ERROR: $BIAS_AA_JSONL not found -- run 15082026_00_pdb_extraction_jsonls_preparation.py first"
    exit 1
fi
if [ ! -f "$PROTEIN_EXPERIMENT_MAP" ]; then
    echo "ERROR: $PROTEIN_EXPERIMENT_MAP not found -- run 15082026_00_pdb_extraction_jsonls_preparation.py first"
    exit 1
fi

mkdir -p "$TABLES_DIR" "$OUTPUTS_DIR"

echo
echo "===== 0. Stage only NOT-yet-designed proteins, stripped to chains A/B ====="
python3 << PYEOF
import csv
import os
import tarfile
from pathlib import Path

STAGE = Path(os.environ["STAGE"])
inputs_prepared = STAGE / "inputs_prepared"
outputs_dir = STAGE / "outputs"
staging_dir = Path(os.environ["MPNN_STAGING_DIR"])
staging_dir.mkdir(parents=True, exist_ok=True)
for child in staging_dir.iterdir():
    if child.is_file():
        child.unlink()

protein_experiment_map = {}
with open(os.environ["PROTEIN_EXPERIMENT_MAP"], newline="") as f:
    for row in csv.DictReader(f):
        protein_experiment_map[row["protein_id"]] = row["experiment"]

GROUP_PREFIXES = {
    'hsl': 'hslarge', 'hsx': 'hsxlarge', 'hss': 'hssmall', 'hsm': 'hsmedium',
    'shl': 'shlarge', 'shm': 'shmedium', 'shx': 'shxlarge', 'shs': 'shsmall',
    'htl': 'htlarge', 'htm': 'htmedium', 'htx': 'htxlarge', 'hts': 'htsmall',
    'thl': 'thlarge', 'thm': 'thmedium', 'thx': 'thxlarge', 'ths': 'thsmall',
}

def group_for_protein(pid):
    prefix = pid[:3]
    group = GROUP_PREFIXES.get(prefix)
    if group is None:
        raise ValueError(f"unrecognized protein_id prefix: {pid!r}")
    return group

def already_designed(pid, experiment, group):
    archive_path = outputs_dir / experiment / f"{group}.tar.gz"
    if not archive_path.exists():
        return False
    with tarfile.open(archive_path, "r:gz") as tf:
        return any(name.startswith(f"{pid}/") for name in tf.getnames())

all_pdbs = sorted(inputs_prepared.glob("*.pdb"))
pending = []
for pdb_path in all_pdbs:
    pid = pdb_path.stem
    group = group_for_protein(pid)
    experiment = protein_experiment_map.get(pid)
    if experiment is None:
        print(f"  [WARNING] {pid}: no experiment mapping in protein_experiment_map.csv "
              f"-- placing under 'unknown_experiment'")
        experiment = "unknown_experiment"
    if already_designed(pid, experiment, group):
        continue
    lines_out = []
    for line in pdb_path.read_text().splitlines():
        if line[:4] in ("ATOM", "HETA"):
            chain = line[21:22].strip()
            if chain not in ("A", "B"):
                continue
        lines_out.append(line)
    (staging_dir / pdb_path.name).write_text("\n".join(lines_out) + "\n")
    pending.append(pid)

print(f"  {len(pending)} protein(s) pending design (of {len(all_pdbs)} prepared total)")
Path(os.environ["PENDING_COUNT_FILE"]).write_text(str(len(pending)))
PYEOF

PENDING_COUNT="$(cat "$PENDING_COUNT_FILE")"
if [ "$PENDING_COUNT" -eq 0 ]; then
    echo
    echo "Nothing new to design -- every prepared protein already has designs in $OUTPUTS_DIR"
    exit 0
fi

echo
echo "===== 1. Parse the complex(es) into ProteinMPNN's jsonl format ====="
python3 "$MPNN_ROOT/helper_scripts/parse_multiple_chains.py" \
    --input_path="$MPNN_STAGING_DIR/" \
    --output_path="$TABLES_DIR/parsed_pdbs.jsonl"

echo
echo "===== 2. Mark BOTH chains as designable (each conditioned on the other) ====="
python3 "$MPNN_ROOT/helper_scripts/assign_fixed_chains.py" \
    --input_path="$TABLES_DIR/parsed_pdbs.jsonl" \
    --output_path="$TABLES_DIR/assigned_pdbs.jsonl" \
    --chain_list "A B"

echo
echo "===== 2b. Tie chain A and chain B together (symmetric design) ====="
python3 "$MPNN_ROOT/helper_scripts/make_tied_positions_dict.py" \
    --input_path="$TABLES_DIR/parsed_pdbs.jsonl" \
    --output_path="$TABLES_DIR/tied_pdbs.jsonl" \
    --homooligomer 1

echo
echo "===== 2c. PRE-FLIGHT: verify symmetric tying on ONE structure before the full batch ====="
verify_symmetry_on_one_structure() {
    local check_dir
    check_dir="$(mktemp -d "$STAGE/.symmetry_check.XXXXXX")"
    trap 'rm -rf "$check_dir"' RETURN
    local test_pdb
    test_pdb="$(ls "$MPNN_STAGING_DIR"/*.pdb 2>/dev/null | head -n 1)"
    if [ -z "$test_pdb" ]; then
        echo "  ERROR: no staged .pdb files found in $MPNN_STAGING_DIR"
        return 1
    fi
    local struct_name
    struct_name="$(basename "$test_pdb" .pdb)"
    echo "  test structure: $struct_name"
    mkdir -p "$check_dir/inputs"
    cp "$test_pdb" "$check_dir/inputs/"
    python3 "$MPNN_ROOT/helper_scripts/parse_multiple_chains.py" \
        --input_path="$check_dir/inputs/" \
        --output_path="$check_dir/parsed.jsonl" >/dev/null
    python3 "$MPNN_ROOT/helper_scripts/assign_fixed_chains.py" \
        --input_path="$check_dir/parsed.jsonl" \
        --output_path="$check_dir/assigned.jsonl" \
        --chain_list "A B" >/dev/null
    python3 "$MPNN_ROOT/helper_scripts/make_tied_positions_dict.py" \
        --input_path="$check_dir/parsed.jsonl" \
        --output_path="$check_dir/tied.jsonl" \
        --homooligomer 1 >/dev/null
    python3 "$MPNN_ROOT/protein_mpnn_run.py" \
        --jsonl_path "$check_dir/parsed.jsonl" \
        --chain_id_jsonl "$check_dir/assigned.jsonl" \
        --fixed_positions_jsonl "$FIXED_POSITIONS_JSONL" \
        --tied_positions_jsonl "$check_dir/tied.jsonl" \
        --bias_AA_jsonl "$BIAS_AA_JSONL" \
        --out_folder "$check_dir/out" \
        --num_seq_per_target 2 \
        --batch_size 2 \
        --sampling_temp "0.2" \
        --backbone_noise 0.1 \
        --seed 37 \
        --path_to_model_weights "$MPNN_ROOT/vanilla_model_weights" \
        --model_name "v_48_020" >/dev/null
    python3 - "$check_dir/out/seqs/${struct_name}.fa" "$struct_name" << 'PYEOF'
import sys
fasta_path, struct_name = sys.argv[1], sys.argv[2]
text = open(fasta_path).read()
records = [r for r in text.split(">") if r.strip()]
designed = records[1:] if len(records) > 1 else records
if not designed:
    print(f"  FAIL: no designed sequences found for {struct_name}")
    sys.exit(1)
all_match = True
for i, record in enumerate(designed):
    header, *seq_lines = record.strip().splitlines()
    seq = "".join(seq_lines)
    chains = seq.split("/")
    if len(chains) != 2:
        print(f"  FAIL: expected 2 chains separated by '/', found {len(chains)} in design {i}")
        all_match = False
        continue
    chain_a, chain_b = chains
    if chain_a == chain_b:
        print(f"  design {i}: chain A == chain B ({len(chain_a)} residues) - OK")
    else:
        all_match = False
        n_diff = sum(1 for x, y in zip(chain_a, chain_b) if x != y)
        print(f"  design {i}: MISMATCH - {n_diff} position(s) differ between chain A and B")
if not all_match:
    print("  RESULT: FAIL - tying is not producing symmetric sequences")
    sys.exit(1)
print("  RESULT: PASS - chain A and B are identical in every design")
PYEOF
    return $?
}
if ! verify_symmetry_on_one_structure; then
    echo
    echo "ERROR: pre-flight symmetry check failed -- aborting before the full batch."
    exit 1
fi
echo "  pre-flight check passed -- continuing to the full batch"

echo
echo "===== 3. Design $NUM_SEQ sequences per structure (tied/symmetric) ====="
python3 "$MPNN_ROOT/protein_mpnn_run.py" \
    --jsonl_path "$TABLES_DIR/parsed_pdbs.jsonl" \
    --chain_id_jsonl "$TABLES_DIR/assigned_pdbs.jsonl" \
    --fixed_positions_jsonl "$FIXED_POSITIONS_JSONL" \
    --tied_positions_jsonl "$TABLES_DIR/tied_pdbs.jsonl" \
    --bias_AA_jsonl "$BIAS_AA_JSONL" \
    --out_folder "$RAW_SEQS_DIR" \
    --num_seq_per_target "$NUM_SEQ" \
    --batch_size "$NUM_SEQ" \
    --sampling_temp "0.2" \
    --backbone_noise 0.1 \
    --seed 37 \
    --path_to_model_weights "$MPNN_ROOT/vanilla_model_weights" \
    --model_name "v_48_020"

echo
echo "===== 4. Split each combined fasta into $NUM_SEQ per-protein files, archive per experiment/group (cumulative) ====="
python3 << PYEOF
import csv
import io
import os
import tarfile
from pathlib import Path

STAGE = Path(os.environ["STAGE"])
NUM_SEQ = $NUM_SEQ
raw_fasta_dir = STAGE / "mpnn_raw" / "seqs"
outputs_dir = STAGE / "outputs"
outputs_dir.mkdir(parents=True, exist_ok=True)

protein_experiment_map = {}
with open(os.environ["PROTEIN_EXPERIMENT_MAP"], newline="") as f:
    for row in csv.DictReader(f):
        protein_experiment_map[row["protein_id"]] = row["experiment"]

GROUP_PREFIXES = {
    'hsl': 'hslarge', 'hsx': 'hsxlarge', 'hss': 'hssmall', 'hsm': 'hsmedium',
    'shl': 'shlarge', 'shm': 'shmedium', 'shx': 'shxlarge', 'shs': 'shsmall',
    'htl': 'htlarge', 'htm': 'htmedium', 'htx': 'htxlarge', 'hts': 'htsmall',
    'thl': 'thlarge', 'thm': 'thmedium', 'thx': 'thxlarge', 'ths': 'thsmall',
}

def group_for_protein(pid):
    return GROUP_PREFIXES[pid[:3]]

def read_archive_members(path):
    members = {}
    if path.exists():
        with tarfile.open(path, "r:gz") as tf:
            for info in tf.getmembers():
                if info.isfile():
                    handle = tf.extractfile(info)
                    members[info.name] = handle.read() if handle else b""
    return members

by_key = {}  # (experiment, group) -> {member_name: bytes}
total_written = 0
n_unknown = 0
for fasta_path in sorted(raw_fasta_dir.glob("*.fa")):
    pid = fasta_path.stem
    group = group_for_protein(pid)
    experiment = protein_experiment_map.get(pid)
    if experiment is None:
        print(f"  [WARNING] {pid}: no experiment mapping in protein_experiment_map.csv "
              f"-- placing under 'unknown_experiment'")
        experiment = "unknown_experiment"
        n_unknown += 1
    text = fasta_path.read_text()
    records = [r for r in text.split(">") if r.strip()]
    designed = records[1:] if len(records) > NUM_SEQ else records
    if len(designed) != NUM_SEQ:
        print(f"  [WARNING] {pid}: expected {NUM_SEQ} designed sequences, found {len(designed)}")
    entries = by_key.setdefault((experiment, group), {})
    for i, record in enumerate(designed[:NUM_SEQ]):
        member_name = f"{pid}/{pid}{i:02d}.fa"
        header, *seq_lines = record.strip().splitlines()
        content = (f">{header}\n" + "\n".join(seq_lines) + "\n").encode()
        entries[member_name] = content
        total_written += 1

for (experiment, group), new_members in sorted(by_key.items()):
    archive_path = outputs_dir / experiment / f"{group}.tar.gz"
    existing = read_archive_members(archive_path)
    collision = set(existing) & set(new_members)
    if collision:
        raise RuntimeError(f"{archive_path}: member name collision on rerun: {sorted(collision)[:5]}")
    merged = {**existing, **new_members}
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = archive_path.parent / f".{archive_path.name}.tmp{os.getpid()}"
    with tarfile.open(tmp_path, "w:gz") as tf:
        for name, content in merged.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            tf.addfile(info, io.BytesIO(content))
    with tarfile.open(tmp_path, "r:gz") as tf:
        names_in_tmp = set(tf.getnames())
    missing = set(merged) - names_in_tmp
    if missing:
        raise RuntimeError(f"{archive_path}: verification failed, missing {sorted(missing)[:5]}")
    os.replace(tmp_path, archive_path)
    print(f"  {experiment}/{group}: +{len(new_members)} file(s) (now {len(merged)} total) -> {archive_path}")

print(f"  wrote {total_written} individual fasta file(s) across {len(by_key)} (experiment, group) archive(s)")
if n_unknown:
    print(f"  [WARNING] {n_unknown} structure(s) had no experiment mapping")
PYEOF

echo
echo "===== 5. Clean up raw ProteinMPNN scratch output (only reached if step 4's archive write + verification above succeeded) ====="
rm -rf "$RAW_SEQS_DIR"
echo "  removed $RAW_SEQS_DIR -- its contents are now safely merged into $OUTPUTS_DIR"

echo
echo "[DONE]"
