#!/usr/bin/env python3
"""
Stage 04 (04_alphafold), script 1: build AlphaFold3 input jsons from the
ProteinMPNN-designed fastas stage 03 already selected and transferred --
inputs/<experiment>/<group>.tar.gz, each holding <protein_id>/<sequence_id>.fa
(exactly 3 per protein, already the best-3 by score -- no re-ranking or
re-selection here, unlike the reference project's 01_af3_json_builder.py,
which pulled straight from a flat outputs/ of ALL designs and had to
re-select the top N itself).

Per sequence, writes exactly 2 AF3 jobs (not 3, unlike the reference
script's a/b/d modes) -- because this project's designs are tied
homooligomers (chain A and chain B are identical by construction, see
15082026_01_run_one_job.sh's make_tied_positions_dict.py --homooligomer 1
step), so a separate "chain B alone" monomer prediction would just
duplicate the "chain A alone" one:
    <sequence_id>_monomer.json  -- single chain "A", the design sequence
    <sequence_id>_dimer.json    -- chains "A" and "B", both tied chains
                                    (whatever the fasta actually recorded
                                    for each, even though they're expected
                                    to match -- see chain-mismatch warning
                                    below)

No ligand/cellulose -- protein-only, same as the reference script.

Output location: 04_alphafold/json/ (flat, not nested under inputs/ or by
experiment) -- sequence_id is globally unique across the whole dataset
(confirmed by every other table in this project keying on it directly),
so no <experiment>/ subfolder is needed to avoid collisions.

Idempotent per job: if both json files for a sequence_id already exist,
that sequence is skipped on rerun -- same "skip what's already done"
convention as every other script here, just checked directly against the
json/ directory's own contents rather than a separate tracking table
(cheap enough that a full existence check every run is fine).

Usage:
    python3 15082026_01_af3_json_builder.py
"""
from __future__ import annotations

import json
import tarfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent          # 04_alphafold/ (when placed in .../scripts/)

INPUTS_DIR = STAGE_ROOT / "inputs"      # <experiment>/<group>.tar.gz, <protein_id>/<sequence_id>.fa (3 per protein)
JSON_DIR = STAGE_ROOT / "json"          # output -- flat, <sequence_id>_monomer.json / <sequence_id>_dimer.json


def discover_group_tars(inputs_dir: Path) -> list[Path]:
    """inputs/<experiment>/<group>.tar.gz -- same <experiment>/<group> layout every stage uses."""
    if not inputs_dir.exists():
        return []
    tars = []
    for experiment_dir in sorted(inputs_dir.iterdir()):
        if not experiment_dir.is_dir():
            continue
        tars.extend(sorted(experiment_dir.glob("*.tar.gz")))
    return tars


def read_member_bytes(archive: tarfile.TarFile, member_name: str) -> bytes:
    handle = archive.extractfile(member_name)
    if handle is None:
        raise FileNotFoundError(f"{member_name} has no extractable content")
    return handle.read()


def split_fasta_chains(content: bytes) -> tuple[str, str]:
    lines = content.decode().strip().splitlines()
    seq_line = lines[1]
    chains = seq_line.split("/")
    if len(chains) != 2:
        raise ValueError(f"expected 2 chains separated by '/', got {len(chains)}: {seq_line!r}")
    return chains[0], chains[1]


def build_af3_json(job_name: str, chains: list[tuple[str, str]]) -> dict:
    """chains: list of (chain_id, sequence). No ligand."""
    return {
        "dialect": "alphafold3",
        "version": 2,
        "name": job_name,
        "sequences": [
            {
                "protein": {
                    "id": chain_id,
                    "sequence": seq,
                    "modifications": [],
                    "unpairedMsa": "",
                    "pairedMsa": "",
                    "templates": [],
                }
            }
            for chain_id, seq in chains
        ],
        "modelSeeds": [1],
        "userCCD": None,
    }


def run_builder() -> None:
    JSON_DIR.mkdir(parents=True, exist_ok=True)

    tar_paths = discover_group_tars(INPUTS_DIR)
    print(f"Building AF3 jsons from {len(tar_paths)} group archive(s) under {INPUTS_DIR}")

    total_sequences = 0
    total_new = 0
    total_skipped = 0
    total_jobs_written = 0
    chain_mismatches = []
    malformed = []

    for tar_path in tar_paths:
        experiment = tar_path.parent.name
        with tarfile.open(tar_path, "r:gz") as archive:
            for name in archive.getnames():
                if "/" not in name or not name.endswith(".fa"):
                    continue
                sequence_id = Path(name).stem  # 'htl000001/htl00000100.fa' -> 'htl00000100'
                total_sequences += 1

                monomer_path = JSON_DIR / f"{sequence_id}_monomer.json"
                dimer_path = JSON_DIR / f"{sequence_id}_dimer.json"
                if monomer_path.exists() and dimer_path.exists():
                    total_skipped += 1
                    continue

                content = read_member_bytes(archive, name)
                try:
                    seq_a, seq_b = split_fasta_chains(content)
                except ValueError as exc:
                    print(f"  [WARNING] {name}: {exc} -- skipped")
                    malformed.append(name)
                    continue
                if seq_a != seq_b:
                    chain_mismatches.append(sequence_id)

                monomer_json = build_af3_json(f"{sequence_id}_monomer", [("A", seq_a)])
                dimer_json = build_af3_json(f"{sequence_id}_dimer", [("A", seq_a), ("B", seq_b)])
                monomer_path.write_text(json.dumps(monomer_json, indent=2))
                dimer_path.write_text(json.dumps(dimer_json, indent=2))
                total_jobs_written += 2
                total_new += 1

    print(f"\nRun complete. Newly built: {total_new} sequence(s) ({total_jobs_written} json file(s)). "
          f"Already-built (skipped): {total_skipped}. Total sequences seen: {total_sequences}.")
    if malformed:
        print(f"[WARNING] {len(malformed)} sequence(s) had an unparseable fasta: {malformed}")
    if chain_mismatches:
        print(f"[WARNING] {len(chain_mismatches)} sequence(s) had chain A != chain B despite the tied-homooligomer "
              f"design -- worth double-checking: {chain_mismatches}")
    print(f"  json dir: {JSON_DIR}")


if __name__ == "__main__":
    run_builder()
