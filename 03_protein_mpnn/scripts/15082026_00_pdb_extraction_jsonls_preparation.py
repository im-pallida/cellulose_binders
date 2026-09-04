#!/usr/bin/env python3
"""
Stage 03, script 0: prepares ProteinMPNN inputs from every group archive
under 03_protein_mpnn/inputs/<experiment>/<group>.tar.gz.

Unlike the old (stage-08-era) version of this script:
  1. Both the .pdb and its companion .json already sit side by side in the
     SAME tar (written by 15082026_05_geometry_filtering.py, carried over
     unchanged by 15082026_02_transfer_to_03.py) -- no separate stage-02
     archive lookup is needed to find a companion file.
  2. The structures are already .pdb (stage 02 writes pdb, not cif), so no
     cif -> pdb conversion is needed either -- the pdb bytes in the tar are
     extracted as-is.
  3. Only the .pdb files are written out (flat, inputs_prepared/), for
     parse_multiple_chains.py to consume. The .json is read in memory only,
     never copied out.
  4. Two ProteinMPNN bias files are written to scripts/jsonls/:
       - fixed_positions.jsonl: per-structure fixed residue positions,
         one combined dict covering every structure prepared so far.
       - bias_AA.jsonl: a single global amino-acid composition bias
         (nudges sampling toward more Ala/Gly), same for every structure.

Fixed positions reuse the exact same fields 15082026_05_geometry_filtering.py
already reads for its own Tyr-anchor alignment (diffused_index_map +
specification.select_exposed): those select_exposed keys ARE the Tyr
positions used as alignment anchors, so they're exactly the residues we
want ProteinMPNN to leave alone. diffused_index_map maps
"<seed chain><seed resi>" -> "<generated chain><generated resi>", and since
stage 02 only rigid-body-aligns (never renumbers), the generated-chain
numbering in that map is already this structure's own final PDB numbering
-- no further re-mapping needed.

Each fixed position is also mirrored onto the OTHER generated chain at the
same local residue number (assign_chains(), copied verbatim from stage 02),
since both physical chains are copies of the same seed monomer sharing one
local numbering scheme (this is the "2xseed_AB" double-ASU design the whole
pipeline is built around).

Incremental/cumulative by design: a protein_id already present in
fixed_positions.jsonl is treated as fully prepared and skipped on rerun,
matching every other script in this pipeline (only ever process what's not
already done; fixed_positions.jsonl doubles as its own tracking table, so
no separate log file is needed).

Also writes scripts/map/protein_experiment_map.csv (protein_id,
experiment) -- CSV, same load-whole/write-whole convention every other
tracking file in this project uses, but kept alongside fixed_positions.jsonl
/ bias_AA.jsonl under scripts/ (its own map/ subfolder) rather than under
tables/, since it's a lookup this stage's own scripts consume, not a
results table. inputs_prepared/ itself is flat (ProteinMPNN's
parse_multiple_chains.py needs one flat input dir), so this table is what
lets 15082026_01_run_*.sh / 15082026_02_3_best_seqs.py sort their own
outputs back into <experiment>/<group>.tar.gz -- without it, the
experiment each protein came from would be lost the moment its pdb lands
in the flat inputs_prepared/ dir.

This table is reconciled on EVERY run, not just when a protein is newly
prepared: every tar is still walked each run regardless of whether its
proteins are already in fixed_positions.jsonl, so a protein_id's
(pid -> experiment) entry gets (re)written every time it's seen, even on
a run where it's otherwise fully skipped. This matters if this table ever
gets out of sync with fixed_positions.jsonl (e.g. deleted/lost separately)
-- rerunning this script re-derives it from the tars without needing to
redo any of the actual pdb-extraction/fixed-position work.
"""
from __future__ import annotations

import csv
import json
import re
import tarfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths -- SCRIPT_DIR/STAGE_ROOT/PROJECT_ROOT convention matches every other
# script in this project, so this stays portable across machines/users.
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent          # 03_protein_mpnn/
PROJECT_ROOT = STAGE_ROOT.parent        # 1cbh_clear/

INPUTS_ROOT = STAGE_ROOT / "inputs"                       # <experiment>/<group>.tar.gz (pdb+json), as written by 15082026_02_transfer_to_03.py
INPUTS_PREPARED_ROOT = STAGE_ROOT / "inputs_prepared"      # flat -- only *.pdb, for parse_multiple_chains.py

JSONLS_DIR = SCRIPT_DIR / "jsonls"
FIXED_POSITIONS_PATH = JSONLS_DIR / "fixed_positions.jsonl"
BIAS_AA_PATH = JSONLS_DIR / "bias_AA.jsonl"

# protein_id -> experiment, so outputs/outputs_clean can be sorted back by
# experiment later even though inputs_prepared/ itself is flat. CSV, kept
# alongside the other scripts/ lookups (jsonls/) rather than in tables/.
MAP_DIR = SCRIPT_DIR / "map"
PROTEIN_EXPERIMENT_MAP_PATH = MAP_DIR / "protein_experiment_map.csv"

# First = seed chain for the "directly named" generated chain (the one
# diffused_index_map actually maps to), second = the mirrored chain --
# same convention/values as 02_geometry_filtering's SEED_CHAINS.
SEED_CHAINS = ('B', 'C')

# 02_geometry_filtering's merge_structures() appends the seed's cellulose
# atoms onto the aligned protein before writing the output pdb (that's why
# real tars carry a 3rd chain like 'M' alongside the 2 protein chains) --
# same chain-letter set cellulose_contacts()/merge_structures() there use.
# Excluded here so assign_chains() below only ever sees the 2 protein
# chains it expects.
CELLULOSE_CHAINS = ('K', 'L', 'M', 'N', 'O')

CHECKPOINT_SIZE = 100

# Global ProteinMPNN amino-acid composition bias (bias_AA.jsonl applies
# uniformly across every position of every structure -- this is NOT a
# per-position bias, that's what fixed_positions.jsonl is for). These are
# placeholder weights -- tune to taste; ProteinMPNN adds this value to that
# amino acid's logit at every design position, so 0.0 = no bias.
BIAS_AA_VALUES = {"A": -1.0, "G": -0.8}


def discover_tars(inputs_root: Path) -> list[Path]:
    """<experiment>/<group>.tar.gz, one level deep, same layout every other stage uses."""
    if not inputs_root.exists():
        return []
    tars = []
    for experiment_dir in sorted(inputs_root.iterdir()):
        if not experiment_dir.is_dir():
            continue
        tars.extend(sorted(experiment_dir.glob("*.tar.gz")))
    return tars


def group_members_by_protein(tar_path: Path) -> dict[str, dict[str, str]]:
    """
    {'htl006254': {'.json': 'htl006254.json', '.pdb': 'htl006254.pdb'}}.
    Filters to .json/.pdb members BEFORE computing protein_id, so a tar's
    own './' directory entry (Path('./').stem == '') never becomes a bogus
    empty-string protein_id -- same fix used in the filtering scripts.
    """
    grouped: dict[str, dict[str, str]] = {}
    with tarfile.open(tar_path, "r:gz") as archive:
        for name in archive.getnames():
            suffix = Path(name).suffix
            if suffix not in (".json", ".pdb"):
                continue
            pid = Path(name).stem
            grouped.setdefault(pid, {})[suffix] = name
    return grouped


def read_member_bytes(archive: tarfile.TarFile, member_name: str) -> bytes:
    handle = archive.extractfile(member_name)
    if handle is None:
        raise FileNotFoundError(f"{member_name} has no extractable content")
    return handle.read()


def pdb_chain_atoms(pdb_bytes: bytes) -> list[dict]:
    """
    Minimal PDB parse -- only pulls each ATOM/HETATM line's chain letter,
    since that's all assign_chains() below needs. Not a general-purpose
    parser (no coordinates, no residue names) on purpose: this script never
    does any geometry, only bookkeeping. Cellulose chains are dropped here
    (see CELLULOSE_CHAINS) so only the 2 protein chains reach assign_chains.
    """
    atoms = []
    for line in pdb_bytes.decode().splitlines():
        if line[:4] in ("ATOM", "HETA"):
            chain = line[21:22].strip()
            if chain in CELLULOSE_CHAINS:
                continue
            atoms.append({'chain': chain})
    return atoms


def split_map_key(key: str) -> tuple[str, int]:
    """Parses a diffused_index_map key/value like 'A24' into ('A', 24). Copied from stage 02."""
    match = re.fullmatch(r"(.+?)(-?\d+)", key)
    if not match:
        raise ValueError(f"cannot parse residue key {key!r}")
    return match.group(1), int(match.group(2))


def assign_chains(atoms: list[dict], diffused_index_map: dict, align_keys: list[str],
                   seed_chains: tuple = SEED_CHAINS) -> dict:
    """
    Copied verbatim (in spirit) from 02_geometry_filtering's assign_chains:
    determines which generated chain is 'directly' named by this protein's
    own align_keys in diffused_index_map, and assigns the other generated
    chain to the other seed chain. Returns {generated_chain: seed_chain}.
    """
    protein_chains = sorted(set(atom['chain'] for atom in atoms))
    if len(protein_chains) != 2:
        raise ValueError(f"expected exactly two generated chains, found {protein_chains}")
    directly_named = set()
    for key in align_keys:
        if key in diffused_index_map:
            chain, _ = split_map_key(diffused_index_map[key])
            directly_named.add(chain)
    if len(directly_named) != 1:
        raise ValueError(
            f"align_keys must all name one directly-mapped generated chain; found {sorted(directly_named)}"
        )
    direct_chain = next(iter(directly_named))
    if direct_chain not in protein_chains:
        raise ValueError(f"json maps chain {direct_chain}, but generated chains are {protein_chains}")
    mirrored_chain = next(c for c in protein_chains if c != direct_chain)
    return {direct_chain: seed_chains[0], mirrored_chain: seed_chains[1]}


def compute_fixed_positions(pdb_bytes: bytes, json_data: dict) -> dict[str, list[int]] | None:
    """
    Returns {'A': [pos, ...], 'B': [pos, ...]} in this structure's own final
    PDB numbering, mirrored across both generated chains, or None if this
    protein's json doesn't carry what's needed (missing map/anchors, or
    chain assignment fails).
    """
    diffused_index_map = json_data.get('diffused_index_map')
    select_exposed = json_data.get('specification', {}).get('select_exposed', '')
    align_keys = [key.strip() for key in select_exposed.split(',') if key.strip()]
    if not diffused_index_map or not align_keys:
        return None
    atoms = pdb_chain_atoms(pdb_bytes)
    chain_assignment = assign_chains(atoms, diffused_index_map, align_keys, SEED_CHAINS)
    direct_chain = next(c for c, s in chain_assignment.items() if s == SEED_CHAINS[0])
    mirrored_chain = next(c for c, s in chain_assignment.items() if s == SEED_CHAINS[1])
    positions = set()
    for key in align_keys:
        if key not in diffused_index_map:
            continue
        _, resi = split_map_key(diffused_index_map[key])
        positions.add(resi)
    if not positions:
        return None
    sorted_positions = sorted(positions)
    return {direct_chain: sorted_positions, mirrored_chain: sorted_positions}


def load_fixed_positions(path: Path) -> dict:
    if not path.exists():
        return {}
    text = path.read_text().strip()
    return json.loads(text) if text else {}


def save_fixed_positions(data: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data) + "\n")


def save_bias_aa(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(BIAS_AA_VALUES) + "\n")


def load_protein_experiment_map(path: Path) -> dict:
    by_id: dict[str, str] = {}
    if path.exists() and path.stat().st_size > 0:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                by_id[row["protein_id"]] = row["experiment"]
    return by_id


def save_protein_experiment_map(data: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["protein_id", "experiment"])
        for pid in sorted(data):
            writer.writerow([pid, data[pid]])


def run_prep() -> None:
    INPUTS_PREPARED_ROOT.mkdir(parents=True, exist_ok=True)
    JSONLS_DIR.mkdir(parents=True, exist_ok=True)
    MAP_DIR.mkdir(parents=True, exist_ok=True)

    fixed_positions_all = load_fixed_positions(FIXED_POSITIONS_PATH)
    already_prepared = len(fixed_positions_all)
    protein_experiment_map = load_protein_experiment_map(PROTEIN_EXPERIMENT_MAP_PATH)

    tar_paths = discover_tars(INPUTS_ROOT)
    print(f"Prep running over {len(tar_paths)} group archive(s) under {INPUTS_ROOT}")
    print(f"Already prepared (in {FIXED_POSITIONS_PATH.name}): {already_prepared}")

    total_new = 0
    total_skipped = 0
    missing_fixed = []
    processed_since_checkpoint = 0

    for tar_path in tar_paths:
        experiment = tar_path.parent.name  # INPUTS_ROOT/<experiment>/<group>.tar.gz
        members = group_members_by_protein(tar_path)
        with tarfile.open(tar_path, "r:gz") as archive:
            for pid, files in sorted(members.items()):
                # Reconciled every run, regardless of already-prepared status
                # below -- this is what lets a lost/out-of-sync
                # protein_experiment_map.csv self-heal on a plain rerun,
                # without needing to redo any pdb-extraction/fixed-position
                # work for proteins that are already done.
                protein_experiment_map[pid] = experiment

                # fixed_positions.jsonl doubles as the "already done" tracker
                # (see module docstring), but it can drift from what's
                # actually sitting in inputs_prepared/ -- e.g. if that
                # directory gets cleaned up/lost separately while
                # fixed_positions.jsonl survives. Requiring the pdb to
                # actually exist on disk too means a pid recorded as done
                # but missing its file here gets re-extracted instead of
                # silently staying missing forever.
                pdb_path = INPUTS_PREPARED_ROOT / f"{pid}.pdb"
                if pid in fixed_positions_all and pdb_path.exists():
                    total_skipped += 1
                    continue
                if ".pdb" not in files:
                    print(f"  [WARNING] {pid}: no .pdb in {tar_path} -- skipping")
                    continue
                if ".json" not in files:
                    print(f"  [WARNING] {pid}: no .json in {tar_path} -- skipping")
                    continue

                pdb_bytes = read_member_bytes(archive, files[".pdb"])
                json_bytes = read_member_bytes(archive, files[".json"])
                pdb_path.write_bytes(pdb_bytes)

                json_data = json.loads(json_bytes.decode())
                try:
                    fixed = compute_fixed_positions(pdb_bytes, json_data)
                except ValueError as exc:
                    print(f"  [WARNING] {pid}: {exc} -- pdb prepared, no fixed positions")
                    fixed = None

                if fixed:
                    fixed_positions_all[pid] = fixed
                else:
                    missing_fixed.append(pid)

                total_new += 1
                processed_since_checkpoint += 1
                if processed_since_checkpoint >= CHECKPOINT_SIZE:
                    save_fixed_positions(fixed_positions_all, FIXED_POSITIONS_PATH)
                    save_protein_experiment_map(protein_experiment_map, PROTEIN_EXPERIMENT_MAP_PATH)
                    print(f"    Checkpoint: {total_new} newly prepared this run, {FIXED_POSITIONS_PATH.name} saved")
                    processed_since_checkpoint = 0

    save_fixed_positions(fixed_positions_all, FIXED_POSITIONS_PATH)
    save_bias_aa(BIAS_AA_PATH)
    save_protein_experiment_map(protein_experiment_map, PROTEIN_EXPERIMENT_MAP_PATH)

    print(f"\nRun complete. Newly prepared this run: {total_new}. Already-prepared (skipped): {total_skipped}.")
    print(f"Total structures with fixed positions recorded: {len(fixed_positions_all)}")
    if missing_fixed:
        print(f"[WARNING] {len(missing_fixed)} newly-prepared structure(s) got a .pdb but NO fixed positions "
              f"(json missing diffused_index_map/select_exposed, or chain assignment failed): {missing_fixed}")
    print(f"  inputs_prepared: {INPUTS_PREPARED_ROOT}")
    print(f"  fixed_positions: {FIXED_POSITIONS_PATH}")
    print(f"  bias_AA:         {BIAS_AA_PATH}")
    print(f"  experiment_map:  {PROTEIN_EXPERIMENT_MAP_PATH}")


if __name__ == "__main__":
    run_prep()
