#!/usr/bin/env python3
"""
Stage 03, script 2: opens every outputs/<experiment>/<group>.tar.gz
archive (written by 15082026_01_run_proteinmpnn.sh -- each holding
<protein_id>/<protein_id>NN.fa, all 16 designed sequences per protein)
and records one row per sequence into tables/stage_03_summary.csv:
    protein_id, experiment, sequence_id, score, ala_content, gly_content,
    other_content

ala_content/gly_content/other_content are the fraction (0.0-1.0) of the
designed sequence that is Ala, Gly, and neither, respectively -- ties back
to step 00's bias_AA.jsonl (which nudges sampling toward more A/G), so
this is how you'd check what composition ProteinMPNN actually produced.
Computed over the full designed sequence with both tied chains combined
(the "/" separator stripped) -- since chain A and chain B are identical
by this pipeline's tied/symmetric design, the fraction is the same
whether computed on one chain or both, so no single-chain special-casing
is needed.

sequence_id is the fasta filename stem (e.g. "htl00000100" for member
htl000001/htl00000100.fa) -- unique on its own across the whole dataset,
so no compound key is needed to track what's already recorded.

score is parsed from the literal "score=" field in each fasta record's
own header, e.g.:
    >T=0.2, sample=7, score=0.51, global_score=0.8734, seq_recovery=0.57
NOTE: this is deliberately the per-residue "score" field, NOT
"global_score" -- the field the rest of this pipeline (best-3 selection,
sorting) has used everywhere else. Confirmed as an explicit choice, not a
default. Any script reading this table downstream should be aware its
"score" column means this field specifically.

Source is outputs/ (every raw design ProteinMPNN produced), not
outputs_clean/ (which only ever holds an already-narrowed top-3) -- this
table is meant to be the complete, unfiltered record that best-3
selection and transfer scripts read from next.

Cumulative by design, same as every other script in this pipeline: a
sequence_id already present in tables/stage_03_summary.csv is treated as
already recorded and skipped on rerun; new rows are merged in and the
table is rewritten in full, sorted by (protein_id, sequence_id).

Usage:
    python3 15082026_02_summary.py
"""
from __future__ import annotations

import csv
import tarfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent          # 03_protein_mpnn/ (when placed in .../scripts/)

OUTPUTS_ROOT = STAGE_ROOT / "outputs"   # <experiment>/<group>.tar.gz, <protein_id>/<protein_id>NN.fa (all 16)
TABLE_PATH = STAGE_ROOT / "tables" / "stage_03_summary.csv"

RESULT_FIELDS = ["protein_id", "experiment", "sequence_id", "score", "ala_content", "gly_content", "other_content"]


def discover_group_tars(outputs_root: Path) -> list[Path]:
    """outputs/<experiment>/<group>.tar.gz -- same <experiment>/<group> layout every other stage uses."""
    if not outputs_root.exists():
        return []
    tars = []
    for experiment_dir in sorted(outputs_root.iterdir()):
        if not experiment_dir.is_dir():
            continue
        tars.extend(sorted(experiment_dir.glob("*.tar.gz")))
    return tars


def group_name_from_tar(tar_path: Path) -> str:
    """Strips the literal '.tar.gz' -- Path.stem only strips ONE suffix ('htlarge.tar'), not both."""
    name = tar_path.name
    if not name.endswith(".tar.gz"):
        raise ValueError(f"expected a .tar.gz file, got {tar_path}")
    return name[: -len(".tar.gz")]


def parse_score(header: str) -> float:
    """
    Header format confirmed from ProteinMPNN's own output, e.g.:
      T=0.2, sample=7, score=0.51, global_score=0.8734, seq_recovery=0.57
    Deliberately matches "score=" only -- "global_score=..." never matches
    since it doesn't start with the literal substring "score=" once split
    on commas and stripped.
    """
    for part in header.split(","):
        part = part.strip()
        if part.startswith("score="):
            return float(part.split("=", 1)[1])
    raise ValueError(f"score not found in header: {header!r}")


def read_member_bytes(archive: tarfile.TarFile, member_name: str) -> bytes:
    handle = archive.extractfile(member_name)
    if handle is None:
        raise FileNotFoundError(f"{member_name} has no extractable content")
    return handle.read()


def compute_composition(seq: str) -> tuple[float, float, float]:
    """
    (ala_content, gly_content, other_content), each a 0.0-1.0 fraction of
    len(seq). seq is expected already-combined across both tied chains
    with the '/' separator stripped (see module docstring).
    """
    if not seq:
        raise ValueError("empty sequence -- cannot compute composition")
    length = len(seq)
    ala_count = seq.count("A")
    gly_count = seq.count("G")
    other_count = length - ala_count - gly_count
    return ala_count / length, gly_count / length, other_count / length


# ---------------------------------------------------------------------------
# Results table -- same load-whole/write-whole CSV-under-tables/ pattern as
# every other tracking table in this project, keyed on sequence_id (unique
# on its own, embeds protein_id as its own prefix).
# ---------------------------------------------------------------------------

def load_table(path: Path) -> dict[str, dict[str, str]]:
    by_sequence_id: dict[str, dict[str, str]] = {}
    if path.exists() and path.stat().st_size > 0:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                by_sequence_id[row["sequence_id"]] = row
    return by_sequence_id


def save_table(by_sequence_id: dict[str, dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_sorted = sorted(by_sequence_id.values(), key=lambda row: (row["protein_id"], row["sequence_id"]))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        for row in rows_sorted:
            writer.writerow({field: row.get(field, "") for field in RESULT_FIELDS})


def run_summary() -> None:
    table = load_table(TABLE_PATH)
    already_done = set(table.keys())

    tar_paths = discover_group_tars(OUTPUTS_ROOT)
    print(f"Summarizing scores across {len(tar_paths)} group archive(s) under {OUTPUTS_ROOT}")
    print(f"Already recorded (in {TABLE_PATH.name}): {len(already_done)}")

    total_new = 0
    total_skipped = 0
    malformed = []

    for tar_path in tar_paths:
        experiment = tar_path.parent.name  # OUTPUTS_ROOT/<experiment>/<group>.tar.gz
        with tarfile.open(tar_path, "r:gz") as archive:
            for name in archive.getnames():
                if "/" not in name or not name.endswith(".fa"):
                    continue
                sequence_id = Path(name).stem  # 'htl000001/htl00000100.fa' -> 'htl00000100'
                if sequence_id in already_done:
                    total_skipped += 1
                    continue
                protein_id = name.split("/", 1)[0]
                content = read_member_bytes(archive, name)
                lines = content.decode().splitlines()
                header = lines[0].lstrip(">")
                seq = "".join(lines[1:]).replace("/", "")
                try:
                    score = parse_score(header)
                    ala_content, gly_content, other_content = compute_composition(seq)
                except ValueError as exc:
                    print(f"  [WARNING] {name}: {exc} -- skipped")
                    malformed.append(name)
                    continue

                table[sequence_id] = {
                    "protein_id": protein_id,
                    "experiment": experiment,
                    "sequence_id": sequence_id,
                    "score": f"{score:.4f}",
                    "ala_content": f"{ala_content:.4f}",
                    "gly_content": f"{gly_content:.4f}",
                    "other_content": f"{other_content:.4f}",
                }
                already_done.add(sequence_id)
                total_new += 1

    if total_new:
        save_table(table, TABLE_PATH)

    print(f"\nRun complete. Newly recorded: {total_new}. Already-recorded (skipped): {total_skipped}.")
    if malformed:
        print(f"[WARNING] {len(malformed)} sequence(s) had an unparseable header: {malformed}")
    print(f"  summary table: {TABLE_PATH} ({len(table)} row(s) total)")


if __name__ == "__main__":
    run_summary()
