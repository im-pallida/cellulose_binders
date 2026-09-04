#!/usr/bin/env python3
"""
Stage 03, script 3 (was 15082026_02_3_best_seqs.py -- renumbered now that
15082026_02_summary.py exists): picks the 3 lowest-"score" sequences per
protein and saves them into outputs_clean/<experiment>/<group>.tar.gz.

Source of truth is now tables/stage_03_summary.csv (written by
15082026_02_summary.py), NOT a fresh header-parse of every fasta record --
this script reads protein_id/experiment/sequence_id/score straight from
that table, groups by protein_id, sorts ascending by score (lower is
better, same convention "score" already uses -- it's ProteinMPNN's own
per-residue negative-log-likelihood, same direction as global_score used
to rank), and takes the best 3.

NOTE on the field switch: earlier versions of this script ranked by
global_score. Per your explicit choice, this recode ranks by "score"
instead -- the same field stage_03_summary.csv records -- so selection no
longer depends on re-parsing fasta headers at all. If you ever want
global_score back as the ranking field, it isn't in stage_03_summary.csv
today and would need to be added there (or re-parsed from the fasta
headers here) -- flag it and I'll wire it up.

The actual .fa file bytes still have to come from the raw
outputs/<experiment>/<group>.tar.gz archive (the summary table only has
metadata, not sequence content) -- <group> is derived from each
protein_id's own 3-character prefix via GROUP_PREFIXES, same scheme every
other per-protein script in this project already uses, so there's no
need to search every tar under an experiment to find the right one.

Cumulative by design, same as every other script in this pipeline: a
protein_id already present in outputs_clean/<experiment>/<group>.tar.gz is
treated as already selected and skipped on rerun; new selections are
merged into the existing group archive (decompress-append-recompress +
atomic replace) rather than overwriting it.

Also writes/updates tables/stage_03_results.csv, one row per protein:
protein_id, experiment, group, time_stamp, and the 3 selected sequences'
names + scores. Rows are written sorted ascending by the best (lowest) of
the 3 scores, same "quality-sorted export" convention as
stage_01_results.csv / stage_02_results.csv.

Usage:
    python3 15082026_03_best_seqs.py
"""
from __future__ import annotations

import csv
import io
import os
import tarfile
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent          # 03_protein_mpnn/ (when placed in .../scripts/)

OUTPUTS_ROOT = STAGE_ROOT / "outputs"              # <experiment>/<group>.tar.gz, <protein_id>/<protein_id>NN.fa (all 16)
OUTPUTS_CLEAN_ROOT = STAGE_ROOT / "outputs_clean"  # <experiment>/<group>.tar.gz, <protein_id>/<protein_id>NN.fa (best 3 only)
SUMMARY_TABLE_PATH = STAGE_ROOT / "tables" / "stage_03_summary.csv"  # input -- written by 15082026_02_summary.py
RESULTS_TABLE_PATH = STAGE_ROOT / "tables" / "stage_03_results.csv"  # output -- this script's own record

BEST_N = 3

# Same 16-entry scheme every other per-protein script in this project
# duplicates locally rather than importing (each stage script is meant to
# run standalone).
GROUP_PREFIXES = {
    'hsl': 'hslarge', 'hsx': 'hsxlarge', 'hss': 'hssmall', 'hsm': 'hsmedium',
    'shl': 'shlarge', 'shm': 'shmedium', 'shx': 'shxlarge', 'shs': 'shsmall',
    'htl': 'htlarge', 'htm': 'htmedium', 'htx': 'htxlarge', 'hts': 'htsmall',
    'thl': 'thlarge', 'thm': 'thmedium', 'thx': 'thxlarge', 'ths': 'thsmall',
}

RESULT_FIELDS = [
    "protein_id",
    "experiment",
    "group",
    "time_stamp",
    "best_seq_1", "best_seq_1_score",
    "best_seq_2", "best_seq_2_score",
    "best_seq_3", "best_seq_3_score",
]


def group_for_protein(protein_id: str) -> str:
    prefix = protein_id[:3]
    group = GROUP_PREFIXES.get(prefix)
    if group is None:
        raise ValueError(f"unrecognized protein_id prefix: {protein_id!r} (from {prefix!r})")
    return group


# ---------------------------------------------------------------------------
# Input: stage_03_summary.csv -- read-only here, grouped by protein_id.
# ---------------------------------------------------------------------------

def load_summary_by_protein(path: Path) -> dict[str, list[dict[str, str]]]:
    """{protein_id: [{'sequence_id':..., 'score':..., 'experiment':...}, ...]}"""
    by_protein: dict[str, list[dict[str, str]]] = {}
    if not path.exists() or path.stat().st_size == 0:
        return by_protein
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            by_protein.setdefault(row["protein_id"], []).append(row)
    return by_protein


# ---------------------------------------------------------------------------
# outputs/ archive access + outputs_clean/ merge -- same pattern as every
# other cumulative archive writer in this project.
# ---------------------------------------------------------------------------

def read_member_bytes(archive: tarfile.TarFile, member_name: str) -> bytes:
    handle = archive.extractfile(member_name)
    if handle is None:
        raise FileNotFoundError(f"{member_name} has no extractable content")
    return handle.read()


def read_archive_members(path: Path) -> dict[str, bytes]:
    members: dict[str, bytes] = {}
    if path.exists():
        with tarfile.open(path, "r:gz") as archive:
            for info in archive.getmembers():
                if info.isfile():
                    handle = archive.extractfile(info)
                    members[info.name] = handle.read() if handle else b""
    return members


def existing_proteins_in_archive(path: Path) -> set[str]:
    """protein_ids already selected in outputs_clean/<experiment>/<group>.tar.gz -- the tracking source, no separate log needed."""
    proteins = set()
    if path.exists():
        with tarfile.open(path, "r:gz") as archive:
            for name in archive.getnames():
                if "/" in name:
                    proteins.add(name.split("/", 1)[0])
    return proteins


def write_merged_archive(archive_path: Path, existing: dict[str, bytes], new_members: dict[str, bytes]) -> None:
    collision = set(existing) & set(new_members)
    if collision:
        raise RuntimeError(f"{archive_path}: member name collision on rerun: {sorted(collision)[:5]}")
    merged = {**existing, **new_members}
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = archive_path.parent / f".{archive_path.name}.tmp{os.getpid()}"
    with tarfile.open(tmp_path, "w:gz") as archive:
        for name, content in merged.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    with tarfile.open(tmp_path, "r:gz") as archive:
        names_in_tmp = set(archive.getnames())
    missing = set(merged) - names_in_tmp
    if missing:
        raise RuntimeError(f"{archive_path}: verification failed, missing {sorted(missing)[:5]}")
    os.replace(tmp_path, archive_path)


# ---------------------------------------------------------------------------
# Results table -- same load-whole/write-whole CSV pattern as every other
# tracking table in this project, keyed on protein_id.
# ---------------------------------------------------------------------------

def load_table(path: Path) -> dict[str, dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    if path.exists() and path.stat().st_size > 0:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                by_id[row["protein_id"]] = row
    return by_id


def save_table(by_id: dict[str, dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # ascending by the best (lowest) of the 3 selected scores -- same
    # "best first" convention as stage_01/02's own sorted_rows.
    rows_sorted = sorted(by_id.values(), key=lambda row: float(row["best_seq_1_score"]))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        for row in rows_sorted:
            writer.writerow({field: row.get(field, "") for field in RESULT_FIELDS})


def build_row(protein_id: str, experiment: str, group: str, time_stamp: str,
              top_entries: list[dict[str, str]]) -> dict[str, str]:
    row = {
        "protein_id": protein_id,
        "experiment": experiment,
        "group": group,
        "time_stamp": time_stamp,
    }
    for rank, entry in enumerate(top_entries, start=1):
        row[f"best_seq_{rank}"] = entry["sequence_id"]
        row[f"best_seq_{rank}_score"] = f"{float(entry['score']):.4f}"
    return row


def run_best_seqs() -> None:
    OUTPUTS_CLEAN_ROOT.mkdir(parents=True, exist_ok=True)

    summary_by_protein = load_summary_by_protein(SUMMARY_TABLE_PATH)
    results_table = load_table(RESULTS_TABLE_PATH)
    print(f"Selecting top {BEST_N} sequences per protein from {SUMMARY_TABLE_PATH} "
          f"({len(summary_by_protein)} protein(s) recorded)")

    total_new_proteins = 0
    total_skipped_proteins = 0
    total_seqs_written = 0
    incomplete_proteins = []
    new_results_rows = 0
    run_time_stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    # Group by (experiment, group) so each source/destination tar is opened once.
    by_archive: dict[tuple[str, str], list[str]] = {}
    for protein_id, rows in summary_by_protein.items():
        experiment = rows[0]["experiment"]
        group = group_for_protein(protein_id)
        by_archive.setdefault((experiment, group), []).append(protein_id)

    for (experiment, group), protein_ids in sorted(by_archive.items()):
        source_path = OUTPUTS_ROOT / experiment / f"{group}.tar.gz"
        clean_archive_path = OUTPUTS_CLEAN_ROOT / experiment / f"{group}.tar.gz"
        if not source_path.exists():
            print(f"  [WARNING] {source_path} not found -- skipping {len(protein_ids)} protein(s) listed for it")
            continue

        already_done = existing_proteins_in_archive(clean_archive_path)
        new_members: dict[str, bytes] = {}

        with tarfile.open(source_path, "r:gz") as archive:
            for protein_id in sorted(protein_ids):
                if protein_id in already_done:
                    total_skipped_proteins += 1
                    continue
                rows = sorted(summary_by_protein[protein_id], key=lambda row: float(row["score"]))
                if len(rows) < BEST_N:
                    incomplete_proteins.append(f"{experiment}/{group}/{protein_id} ({len(rows)}/{BEST_N} scored)")
                top_rows = rows[:BEST_N]
                top_entries = []
                for row in top_rows:
                    member_name = f"{protein_id}/{row['sequence_id']}.fa"
                    content = read_member_bytes(archive, member_name)
                    new_members[member_name] = content
                    top_entries.append(row)
                if top_entries:
                    results_table[protein_id] = build_row(protein_id, experiment, group, run_time_stamp, top_entries)
                    new_results_rows += 1
                total_new_proteins += 1

        if not new_members:
            continue
        existing = read_archive_members(clean_archive_path)
        write_merged_archive(clean_archive_path, existing, new_members)
        total_seqs_written += len(new_members)
        print(f"  {experiment}/{group}: {len(new_members)} sequence(s) across newly-selected protein(s) "
              f"(now {len(existing) + len(new_members)} total) -> {clean_archive_path}")

    if new_results_rows:
        save_table(results_table, RESULTS_TABLE_PATH)

    print(f"\nRun complete. Newly selected: {total_new_proteins} protein(s) ({total_seqs_written} sequence(s) written). "
          f"Already-selected (skipped): {total_skipped_proteins}.")
    if incomplete_proteins:
        print(f"[WARNING] {len(incomplete_proteins)} protein(s) had fewer than {BEST_N} scored sequences: "
              f"{incomplete_proteins}")
    print(f"  outputs_clean:  {OUTPUTS_CLEAN_ROOT}")
    print(f"  results table:  {RESULTS_TABLE_PATH} ({len(results_table)} row(s) total)")


if __name__ == "__main__":
    run_best_seqs()
