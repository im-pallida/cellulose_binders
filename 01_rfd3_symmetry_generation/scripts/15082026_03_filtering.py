#!/usr/bin/env python3
"""
Stage 01 filter: n_chainbreaks + non_loop_fraction (radius_of_gyration: tracked only)

Combines two things you already had, into one script:
  - the biological criteria + per-metric error handling from
    22072026_03_filter_rfd3_metrics_one.py (stage02)
  - the tar-scanning / passed-rejected-file-routing / sorted-table
    mechanism from 31072026_geometry_filtering.py

I/O
---
input:
  -- <stage_root>/outputs_clean/<experiment>/*.tar.gz
     each tar contains, per protein_id: "<protein_id>.json" (metadata
     with a "metrics" object) and "<protein_id>.cif" or "<protein_id>.pdb"
     (structure file)
output:
  -- <stage_root>/tables/stage_01_results.csv
     one row per protein_id, PASSED / REJECTED / ERROR
  -- <stage_root>/outputs_clean_stage01/<experiment>/{passed,rejected}/<experiment>_<outcome>.tar.gz
     the json+structure pair for every PASSED/REJECTED protein, regrouped
  -- printed per-protein output as each tar is scanned, and a running
     counts summary at the end

Rule (biological, same cutoffs as stage02)
-------------------------------------------
  -- n_chainbreaks == 0                    -> else REJECT (chainbreaks)
  -- non_loop_fraction > NON_LOOP_CUTOFF   -> else REJECT (loopy_structure)
  -- radius_of_gyration: recorded in the table only, never rejects

ERROR vs REJECTED (same distinction as stage02)
------------------------------------------------
  REJECTED = valid metrics, failed the biological cutoffs.
  ERROR    = the protein couldn't be evaluated at all (missing json,
             missing structure file, malformed json, missing/non-numeric
             metric). Recorded in the table with an error_reason so it's
             not silently dropped, but its files are left untouched in
             the source tar (nothing to safely route) and it is skipped
             on future runs the same as a PASSED/REJECTED row would be.

Already-processed proteins (protein_id already present in
stage_01_results.csv) are skipped without re-extracting them from the
tar -- same idea as geometry_filtering.py's load_table/skip check, but
keyed on protein_id alone since this filter has no cutoff sweep.

Known limitation (inherited from geometry_filtering.py's own
regroup_and_archive): each run writes ONE fresh
<experiment>_<outcome>.tar.gz from that run's newly-routed loose files.
If you re-run after manually clearing rows out of the table, the new
archive will NOT contain protein files archived by a previous run --
only this run's. If you need cross-run merging, use
pipeline_tracking.append_to_growing_archive's decompress-append-recompress
pattern instead of the plain tarfile.open(path, "w:gz") below.

Layers
------
1. I/O: archive_paths, group_members_by_protein, extract_member,
   route_files, regroup_and_archive, load_table, save_table
2. Metrics + criteria: numeric_metric, read_metrics, decide_status
3. Row assembly: build_row
4. Orchestrator: run_filter / main
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
import os
import sys
import tarfile
from pathlib import Path
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# CONSTANTS -- adjust to match your real directory tree before running.
# Assumes this script lives in <stage_root>/scripts/, same convention as
# geometry_filtering.py's STAGE04 = Path(__file__).resolve().parent.parent
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent # scripts dir
STAGE_ROOT = SCRIPT_DIR.parent # stage dir

INPUTS_ROOT = STAGE_ROOT / "outputs_clean"                    # tar.gz archives to scan, one subdir per experiment
OUTPUTS_RAW_ROOT = STAGE_ROOT / "sorted_raw"          # scratch: loose passed/rejected files before re-archiving
OUTPUTS_CLEAN_ROOT = STAGE_ROOT / "sorted_clean"      # final passed/rejected tar.gz per experiment
TABLES_DIR = STAGE_ROOT / "tables"
TABLE_PATH = TABLES_DIR / "stage_01_results.csv"

NON_LOOP_CUTOFF = 0.7
CHECKPOINT_SIZE = 100

RESULT_FIELDS = [
    "protein_id",
    "experiment_name",
    "time_stamp",
    "n_chainbreaks",
    "non_loop_fraction",
    "non_loop_cutoff",
    "radius_of_gyration",
    "status",
    "rejection_reason",
    "error_reason",
]

# Same scheme as geometry_filtering.py's GROUP_PREFIXES: first 3 characters
# of protein_id -> full group name used for archive naming.
GROUP_PREFIXES = {
    'hsl': 'hslarge',
    'hsx': 'hsxlarge',
    'hss': 'hssmall',
    'hsm': 'hsmedium',
    'shl': 'shlarge',
    'shm': 'shmedium',
    'shx': 'shxlarge',
    'shs': 'shsmall',
    'htl': 'htlarge',
    'htm': 'htmedium',
    'htx': 'htxlarge',
    'hts': 'htsmall',
    'thl': 'thlarge',
    'thm': 'thmedium',
    'thx': 'thxlarge',
    'ths': 'thsmall',
}

class Stage01InputError(Exception):
    """Malformed/missing input or metric -- pipeline error, not a biological rejection."""


# ---------------------------------------------------------------------------
# Layer 1: I/O
# ---------------------------------------------------------------------------

def collect_dirs(root_dir: Path) -> list[Path]: # in the future I will use pathlib instead of os everywhere; func is aiming to collect the dirs
    return sorted(p for p in root_dir.iterdir() if p.is_dir()) # good check for dirs


def archive_paths(inputs_root: Path, specific_experiment: str | None = None, # func for filtering only the specific tar/directory
                   specific_tar_gz: str | None = None) -> list[Path]:
    """Same shape as stage 02 archive_direction, Path-based."""
    if specific_tar_gz is not None:
        return [inputs_root / specific_tar_gz]

    if specific_experiment is not None:
        experiments = [inputs_root / specific_experiment]
    else:
        experiments = collect_dirs(inputs_root)

    tars: list[Path] = []
    for experiment in experiments:
        if experiment.is_dir():
            tars.extend(sorted(experiment.glob("*.tar.gz"))) # glob helps to find matching patterns
    return tars


def protein_id_from_member(member_name: str) -> str:
    return Path(member_name).stem # extracts the file name without the extension


def experiment_name_from_tar(tar_path: Path) -> str:
    return tar_path.parent.name


def group_members_by_protein(member_names: list[str]) -> dict[str, dict[str, str]]:
    """
    Groups tar members by protein_id -> {'json': name, 'structure': name}.
    A member with an unmatched suffix (neither .json/.cif/.pdb) is
    skipped entirely -- important because tar archives built with e.g.
    `tar czf out.tar.gz -C dir .` include a directory entry ("./") whose
    Path(...).stem is "", which would otherwise create a bogus
    empty-protein_id group. Either of json/structure can still be
    missing for a real protein_id -- that's decided/handled in
    run_filter, not here, so grouping itself never raises.
    """
    groups: dict[str, dict[str, str]] = {}
    for member_name in member_names:
        suffix = Path(member_name).suffix.lower()
        if suffix not in (".json", ".cif", ".pdb"):
            continue
        pid = protein_id_from_member(member_name)
        entry = groups.setdefault(pid, {})
        if suffix == ".json":
            entry["json"] = member_name
        else:
            entry["structure"] = member_name
    return groups


def extract_member(archive: tarfile.TarFile, member_name: str, extract_dir: Path) -> Path:
    extract_dir.mkdir(parents=True, exist_ok=True)
    archive.extract(member_name, path=extract_dir)
    return extract_dir / member_name


def route_files(json_path: Path, structure_path: Path, experiment_raw_dir: Path,
                 outcome_dir_name: str) -> None:
    """
    Moves the already-extracted json + structure file into
    outputs_raw_stage01/<experiment>/<passed|rejected>/ -- loose files,
    grouped into one tar.gz later by regroup_and_archive. Same
    loose-then-archive two-step as geometry_filtering.py, so a crash
    mid-run never leaves a half-written tar.gz.
    """
    target_dir = experiment_raw_dir / outcome_dir_name
    target_dir.mkdir(parents=True, exist_ok=True)
    json_path.rename(target_dir / json_path.name)
    structure_path.rename(target_dir / structure_path.name)


def group_for_protein(protein_id: str) -> str:
    """Same scheme as geometry_filtering.py's group_for_protein: the first
    3 characters of protein_id encode order+size (e.g. 'htm' -> 'htmedium')."""
    prefix = protein_id[:3]
    if prefix not in GROUP_PREFIXES:
        raise ValueError(
            f"Unknown group prefix {prefix!r} for protein_id {protein_id!r}. "
            f"Expected one of: {sorted(GROUP_PREFIXES)}"
        )
    return GROUP_PREFIXES[prefix]


def read_archive_members(archive_path: Path) -> list[tuple[tarfile.TarInfo, bytes]]:
    members: list[tuple[tarfile.TarInfo, bytes]] = []
    with tarfile.open(archive_path, "r:gz") as archive:
        for member_info in archive.getmembers():
            if not member_info.isfile():
                continue
            extracted = archive.extractfile(member_info)
            if extracted is None:
                raise RuntimeError(f"could not read member {member_info.name!r} from {archive_path}")
            members.append((member_info, extracted.read()))
    return members


def regroup_and_archive(experiment_raw_dir: Path, experiment_clean_dir: Path) -> None:
    """
    Buckets each outcome's loose files by group (group_for_protein) and
    MERGES them into <group>.tar.gz per bucket under
    experiment_clean_dir/<outcome>/ -- matches outputs_clean's own
    per-group naming (htmedium.tar.gz, thlarge.tar.gz, ...) instead of
    one flat <experiment>_<outcome>.tar.gz.
    Cumulative across runs: any members already in an existing
    <group>.tar.gz are read first and kept, this run's new loose files
    are added, the combined result is written to a temp file, verified
    complete, and atomically swapped in (os.replace) -- only then are
    the loose originals deleted. Same decompress-append-recompress
    approach the sync scripts (15082026_04_transfer_to_02.py) already
    use on the destination side, applied here at the source so
    sorted_clean/ itself is a complete history, not just this run's
    batch. Scoped to one group at a time so a bad protein_id prefix in
    one group doesn't stop the other groups in the same outcome from
    archiving successfully.
    """
    for outcome in ("passed", "rejected"):
        raw_outcome_dir = experiment_raw_dir / outcome
        if not raw_outcome_dir.is_dir():
            continue
        loose_files = sorted(p for p in raw_outcome_dir.iterdir() if p.is_file())
        if not loose_files:
            continue

        files_by_group: dict[str, list[Path]] = {}
        for file_path in loose_files:
            group = group_for_protein(file_path.stem)
            files_by_group.setdefault(group, []).append(file_path)

        clean_outcome_dir = experiment_clean_dir / outcome
        clean_outcome_dir.mkdir(parents=True, exist_ok=True)

        for group, group_files in sorted(files_by_group.items()):
            archive_path = clean_outcome_dir / f"{group}.tar.gz"

            existing_members: list[tuple[tarfile.TarInfo, bytes]] = []
            if archive_path.exists():
                existing_members = read_archive_members(archive_path)
            existing_names = {member_info.name for member_info, _ in existing_members}

            new_names = {file_path.name for file_path in group_files}
            collision = existing_names & new_names
            if collision:
                raise RuntimeError(
                    f"{archive_path}: these members are already in the archive but were "
                    f"about to be re-added: {sorted(collision)}"
                )

            temp_path = clean_outcome_dir / f".{group}.tar.gz.tmp{os.getpid()}"
            with tarfile.open(temp_path, "w:gz") as archive:
                for member_info, payload in existing_members:
                    archive.addfile(member_info, io.BytesIO(payload))
                for file_path in group_files:
                    archive.add(file_path, arcname=file_path.name)

            with tarfile.open(temp_path, "r:gz") as archive:
                archived_names = {Path(name).name for name in archive.getnames()}

            expected_names = existing_names | new_names
            missing = expected_names - archived_names
            if missing:
                temp_path.unlink()
                raise RuntimeError(f"Verification failed for {archive_path}: missing {sorted(missing)}")

            os.replace(temp_path, archive_path)

            for file_path in group_files:
                file_path.unlink()


def load_table(table_path: Path) -> dict[str, dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    if table_path.exists() and table_path.stat().st_size > 0:
        with table_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                by_id[row["protein_id"]] = row
    return by_id


def save_table(by_id: dict[str, dict[str, str]], table_path: Path) -> None:
    table_path.parent.mkdir(parents=True, exist_ok=True)
    # worst non_loop_fraction first (blank/error rows sort first via "" < any float string
    # only by luck of formatting, so sort on a parsed key instead)
    def sort_key(row: dict[str, str]) -> tuple:
        try:
            nlf = float(row.get("non_loop_fraction", ""))
        except ValueError:
            nlf = float("-inf")  # ERROR rows (no numeric value) sort first, for visibility
        return (row.get("status") != "ERROR", nlf, row.get("protein_id", ""))

    rows_sorted = sorted(by_id.values(), key=sort_key)
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        for row in rows_sorted:
            writer.writerow({field: row.get(field, "") for field in RESULT_FIELDS})


def print_counts(by_id: dict[str, dict[str, str]]) -> None:
    print()
    print("===== STAGE 01 RUNNING COUNTS =====")
    print(f"total tracked = {len(by_id)}")
    counts: dict[str, int] = {}
    for row in by_id.values():
        status = row.get("status", "")
        if status:
            counts[status] = counts.get(status, 0) + 1
    if counts:
        print("  " + ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))


# ---------------------------------------------------------------------------
# Layer 2: metrics + criteria
# ---------------------------------------------------------------------------

def numeric_metric(metrics: dict[str, Any], name: str) -> float:
    if name not in metrics:
        raise Stage01InputError(f"missing_metric:{name}")
    value = metrics[name]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise Stage01InputError(f"non_numeric_metric:{name}")
    value = float(value)
    if not math.isfinite(value):
        raise Stage01InputError(f"non_finite_metric:{name}")
    return value


def read_metrics(json_path: Path) -> dict[str, float]:
    try:
        with json_path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
    except json.JSONDecodeError as exc:
        raise Stage01InputError(f"malformed_json:line_{exc.lineno}:column_{exc.colno}") from exc
    except OSError as exc:
        raise Stage01InputError(f"unreadable_json:{exc}") from exc

    metrics = payload.get("metrics") if isinstance(payload, dict) else None
    if not isinstance(metrics, dict):
        raise Stage01InputError("missing_metrics_object")

    return {
        "n_chainbreaks": numeric_metric(metrics, "n_chainbreaks"),
        "non_loop_fraction": numeric_metric(metrics, "non_loop_fraction"),
        "radius_of_gyration": numeric_metric(metrics, "radius_of_gyration"),
    }


def decide_status(n_chainbreaks: float, non_loop_fraction: float) -> tuple[str, list[str]]:
    reasons = []
    if n_chainbreaks != 0:
        reasons.append("chainbreaks")
    if not (non_loop_fraction > NON_LOOP_CUTOFF):
        reasons.append("loopy_structure")
    status = "REJECTED" if reasons else "PASSED"
    return status, reasons


# ---------------------------------------------------------------------------
# Layer 3: row assembly
# ---------------------------------------------------------------------------

def build_row(protein_id: str, experiment_name: str, metrics: dict[str, float],
              status: str, reasons: list[str], error_reason: str = "") -> dict[str, str]:
    return {
        "protein_id": protein_id,
        "experiment_name": experiment_name,
        "time_stamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_chainbreaks": f'{metrics["n_chainbreaks"]:.6g}' if "n_chainbreaks" in metrics else "",
        "non_loop_fraction": f'{metrics["non_loop_fraction"]:.10g}' if "non_loop_fraction" in metrics else "",
        "non_loop_cutoff": f"{NON_LOOP_CUTOFF:.10g}",
        "radius_of_gyration": f'{metrics["radius_of_gyration"]:.10g}' if "radius_of_gyration" in metrics else "",
        "status": status,
        "rejection_reason": ";".join(reasons),
        "error_reason": error_reason,
    }


# ---------------------------------------------------------------------------
# Layer 4: orchestrator
# ---------------------------------------------------------------------------

def run_filter(specific_experiment: str | None = None, specific_tar_gz: str | None = None) -> int:
    if not INPUTS_ROOT.is_dir():
        print(f"RESULT: ERROR (inputs_root not found: {INPUTS_ROOT.resolve()})")
        return 1

    by_id = load_table(TABLE_PATH)

    tar_paths = archive_paths(INPUTS_ROOT, specific_experiment, specific_tar_gz)

    print("===== STAGE 01 FILTER: chainbreaks + non_loop_fraction =====")
    print(f"inputs_root  = {INPUTS_ROOT.resolve()}")
    print(f"table_path   = {TABLE_PATH.resolve()}")
    print(f"criteria     = n_chainbreaks == 0, non_loop_fraction > {NON_LOOP_CUTOFF} "
          f"(radius_of_gyration: recorded only)")
    print()

    touched_experiments: dict[str, Path] = {}
    processed_since_checkpoint = 0

    for tar_path in tar_paths:
        experiment_name = experiment_name_from_tar(tar_path)
        experiment_raw_dir = OUTPUTS_RAW_ROOT / experiment_name
        scratch_dir = experiment_raw_dir / "_scratch"

        with tarfile.open(tar_path) as archive:
            groups = group_members_by_protein(archive.getnames())

            for protein_id, members in groups.items():
                if protein_id in by_id:
                    print(f"{protein_id} skipped (already in table)")
                    continue

                json_member = members.get("json")
                structure_member = members.get("structure")

                if json_member is None:
                    by_id[protein_id] = build_row(protein_id, experiment_name, {}, "ERROR", [],
                                                   error_reason="missing_metadata_json")
                    print(f"{protein_id}: ERROR (missing_metadata_json)")
                    continue
                if structure_member is None:
                    by_id[protein_id] = build_row(protein_id, experiment_name, {}, "ERROR", [],
                                                   error_reason="missing_structure_file")
                    print(f"{protein_id}: ERROR (missing_structure_file)")
                    continue

                json_path = extract_member(archive, json_member, scratch_dir)
                structure_path = extract_member(archive, structure_member, scratch_dir)

                try:
                    metrics = read_metrics(json_path)
                except Stage01InputError as exc:
                    by_id[protein_id] = build_row(protein_id, experiment_name, {}, "ERROR", [],
                                                   error_reason=str(exc))
                    print(f"{protein_id}: ERROR ({exc})")
                    json_path.unlink(missing_ok=True)
                    structure_path.unlink(missing_ok=True)
                    continue

                status, reasons = decide_status(metrics["n_chainbreaks"], metrics["non_loop_fraction"])
                by_id[protein_id] = build_row(protein_id, experiment_name, metrics, status, reasons)

                print(f"{protein_id}: n_chainbreaks={metrics['n_chainbreaks']:.4g} "
                      f"non_loop_fraction={metrics['non_loop_fraction']:.4g} "
                      f"radius_of_gyration={metrics['radius_of_gyration']:.4g} (tracked only) "
                      f"-> {status}" + (f" ({';'.join(reasons)})" if reasons else ""))

                outcome_dir_name = "passed" if status == "PASSED" else "rejected"
                route_files(json_path, structure_path, experiment_raw_dir, outcome_dir_name)
                touched_experiments[experiment_name] = experiment_raw_dir

                processed_since_checkpoint += 1
                if processed_since_checkpoint >= CHECKPOINT_SIZE:
                    save_table(by_id, TABLE_PATH)
                    print(f"Checkpoint: {processed_since_checkpoint} proteins processed this run, table saved.")
                    processed_since_checkpoint = 0

    save_table(by_id, TABLE_PATH)

    for experiment_name, experiment_raw_dir in touched_experiments.items():
        regroup_and_archive(experiment_raw_dir, OUTPUTS_CLEAN_ROOT / experiment_name)

    print_counts(by_id)
    print(f"Run complete. Total proteins in table: {len(by_id)}")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--experiment", default=None,
                         help="Only scan this one experiment subfolder of outputs_clean/")
    parser.add_argument("--tar-gz", default=None,
                         help="Only scan this one tar.gz (path relative to outputs_clean/)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    return run_filter(specific_experiment=args.experiment, specific_tar_gz=args.tar_gz)


if __name__ == "__main__":
    sys.exit(main())
