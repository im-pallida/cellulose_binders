#!/usr/bin/env python3
"""
Stage 04 (04_alphafold), script 3: copies the ProteinMPNN INPUT backbone
(03_protein_mpnn/inputs_prepared/<protein_id>.pdb -- the exact structure
each design's sequence was generated to fit) into
04_alphafold/inputs/reference/<protein_id>.pdb, for every protein_id that
actually has sequences transferred into this stage (inputs/<experiment>/
<group>.tar.gz). 15082026_04_compute_rmsd.py uses these as the ground
truth each AF3 prediction is compared against.

Only pulls references for protein_ids we actually need (discovered from
this stage's own inputs/ tars), not the whole of 03_protein_mpnn's
inputs_prepared/ -- keeps this directory scoped to what stage 04 is
actually working on.

Cumulative by design, same as every other transfer script in this
project: a protein_id whose reference pdb is already present on disk is
skipped on rerun.

Usage:
    python3 15082026_03_transfer_reference.py
"""
from __future__ import annotations

import shutil
import tarfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent               # 04_alphafold/
PROJECT_ROOT = STAGE_ROOT.parent

INPUTS_DIR = STAGE_ROOT / "inputs"           # <experiment>/<group>.tar.gz, <protein_id>/<sequence_id>.fa
REFERENCE_DIR = INPUTS_DIR / "reference"     # output -- flat, <protein_id>.pdb
SOURCE_DIR = PROJECT_ROOT / "03_protein_mpnn" / "inputs_prepared"  # flat, <protein_id>.pdb


def discover_group_tars(inputs_dir: Path) -> list[Path]:
    """inputs/<experiment>/<group>.tar.gz -- same <experiment>/<group> layout every stage uses.
    Deliberately does NOT descend into inputs/reference/ -- that's this
    script's own output, not a group archive."""
    if not inputs_dir.exists():
        return []
    tars = []
    for experiment_dir in sorted(inputs_dir.iterdir()):
        if not experiment_dir.is_dir() or experiment_dir.name == "reference":
            continue
        tars.extend(sorted(experiment_dir.glob("*.tar.gz")))
    return tars


def needed_protein_ids(tar_paths: list[Path]) -> set[str]:
    needed = set()
    for tar_path in tar_paths:
        with tarfile.open(tar_path, "r:gz") as archive:
            for name in archive.getnames():
                if "/" not in name:
                    continue
                needed.add(name.split("/", 1)[0])
    return needed


def run_transfer() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    tar_paths = discover_group_tars(INPUTS_DIR)
    needed = needed_protein_ids(tar_paths)
    print(f"Found {len(needed)} protein_id(s) needing a reference pdb, across {len(tar_paths)} group archive(s)")

    total_new = 0
    total_skipped = 0
    missing = []

    for protein_id in sorted(needed):
        dest_path = REFERENCE_DIR / f"{protein_id}.pdb"
        if dest_path.exists():
            total_skipped += 1
            continue
        source_path = SOURCE_DIR / f"{protein_id}.pdb"
        if not source_path.exists():
            print(f"  [WARNING] {protein_id}: no reference pdb found at {source_path} -- skipping")
            missing.append(protein_id)
            continue
        shutil.copy2(source_path, dest_path)
        total_new += 1

    print(f"\nRun complete. Newly copied: {total_new}. Already present (skipped): {total_skipped}.")
    if missing:
        print(f"[WARNING] {len(missing)} protein_id(s) had no reference pdb in {SOURCE_DIR}: {missing}")
    print(f"  reference dir: {REFERENCE_DIR}")


if __name__ == "__main__":
    run_transfer()
