#!/usr/bin/env python3
"""
Stage 04 (04_alphafold), script 4: for every sequence with BOTH its
monomer and dimer AF3 jobs done, computes a ChimeraX-matchmaker-style
RMSD (sequence-align to establish residue correspondence, CA-atom Kabsch
superposition, iterative outlier pruning) against the ProteinMPNN INPUT
backbone that sequence was designed to fit (inputs/reference/<protein_id>.pdb,
written by 15082026_03_transfer_reference.py).

Monomer job: compares the single predicted chain "A" against the
reference's chain A.
Dimer job: does ONE JOINT superposition across BOTH chains' matched CA
atoms together, so the result reflects whether the predicted RELATIVE
chain arrangement (association) matches the reference, not just each
chain's internal fold accuracy.

Of AF3's 5 predicted samples per job, only the one AF3 itself ranks #1
(via ranking_scores.csv) is scored -- that's the prediction you'd
actually use downstream, so pass/fail reflects that one, not whichever of
the 5 happens to look best in isolation.

Pass/fail: RMSD_THRESHOLD_MONOMER/DIMER and
MATCHED_FRACTION_THRESHOLD_MONOMER/DIMER below are separate constants for
monomer vs dimer (both 5.0 / 0.8 for now, per your instruction) --
matched_fraction is matched_pairs/total_aligned_pairs, i.e. "matched
pairs > 80%". A sequence PASSES only if BOTH its monomer and dimer jobs
individually pass.

Writes/updates tables/stage_04_results.csv -- cumulative, same convention
as every other results table in this project: a sequence_id already
present is skipped on rerun; a sequence whose monomer+dimer AF3 jobs
aren't both finished yet, or whose reference pdb hasn't been transferred
yet, is left for a future run (not recorded as a failure).

Also sorts each newly-scored sequence's full AF3 output (both monomer and
dimer job folders -- model.cif, ranking_scores.csv, everything AF3wrote)
into outputs_clean/<passed|rejected>/<experiment>/<group>.tar.gz, nested
under <protein_id>/<sequence_id>_monomer/... and
<protein_id>/<sequence_id>_dimer/... -- cumulative merge into the
existing archive, same decompress-append-recompress + atomic-replace
pattern as every other archive writer here.

AF3's exact output directory/file naming is assumed to match your
confirmed-working reference project (<job_dir>/seed-*_sample-*/model.cif
+ ranking_scores.csv) -- this script searches defensively and prints
exactly what it found, so any mismatch is visible immediately.

Usage:
    python3 15082026_04_compute_rmsd.py
"""
from __future__ import annotations

import csv
import io
import os
import tarfile
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent               # 04_alphafold/

INPUTS_DIR = STAGE_ROOT / "inputs"           # <experiment>/<group>.tar.gz, <protein_id>/<sequence_id>.fa
REFERENCE_DIR = INPUTS_DIR / "reference"     # <protein_id>.pdb, written by 15082026_03_transfer_reference.py
OUTPUTS_DIR = STAGE_ROOT / "outputs"         # <sequence_id>_monomer/, <sequence_id>_dimer/ (AF3 job output)
OUTPUTS_CLEAN_DIR = STAGE_ROOT / "outputs_clean"  # <passed|rejected>/<experiment>/<group>.tar.gz
RESULTS_TABLE_PATH = STAGE_ROOT / "tables" / "stage_04_results.csv"

# ---------------------------------------------------------------------
# Pass/fail thresholds -- separate for monomer vs dimer per your
# instruction, both 5.0 / 0.8 for now.
# ---------------------------------------------------------------------
RMSD_THRESHOLD_MONOMER = 5.0
MATCHED_FRACTION_THRESHOLD_MONOMER = 0.8
RMSD_THRESHOLD_DIMER = 5.0
MATCHED_FRACTION_THRESHOLD_DIMER = 0.8

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

RESULT_FIELDS = [
    "protein_id", "experiment", "group", "sequence_id",
    "monomer_rmsd", "monomer_matched_pairs", "monomer_total_aligned_pairs",
    "monomer_matched_fraction", "monomer_ranking_score", "monomer_pass",
    "dimer_rmsd", "dimer_matched_pairs", "dimer_total_aligned_pairs",
    "dimer_matched_fraction", "dimer_ranking_score", "dimer_pass",
    "passed", "destination",
]


# ---------------------------------------------------------------------
# mmCIF parsing (generic -- reads the actual column header, doesn't
# assume a fixed column order)
# ---------------------------------------------------------------------
def parse_mmcif_atom_site(cif_text):
    lines = cif_text.splitlines()
    col_names = []
    data_start = None
    for i, line in enumerate(lines):
        if line.strip() == "loop_":
            j = i + 1
            peek_cols = []
            while j < len(lines) and lines[j].strip().startswith("_atom_site."):
                peek_cols.append(lines[j].strip()[len("_atom_site."):])
                j += 1
            if peek_cols:
                col_names = peek_cols
                data_start = j
                break
    if not col_names:
        raise ValueError("no _atom_site loop_ found in this mmCIF text")
    name_to_idx = {name: idx for idx, name in enumerate(col_names)}
    rows = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("_") or stripped == "loop_":
            break
        fields = stripped.split()
        if fields[0] not in ("ATOM", "HETATM"):
            break
        rows.append({name: fields[name_to_idx[name]] for name in col_names if name_to_idx[name] < len(fields)})
    return rows


def extract_chain_ca_mmcif(rows, chain_id):
    by_resnum = {}
    for r in rows:
        if r.get("auth_asym_id") != chain_id:
            continue
        if r.get("label_atom_id") != "CA":
            continue
        resnum = int(r["auth_seq_id"])
        by_resnum[resnum] = (r.get("label_comp_id"), (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])))
    ordered = sorted(by_resnum.items())
    seq = "".join(AA3TO1.get(resname, "X") for _, (resname, _) in ordered)
    coords = np.array([xyz for _, (_, xyz) in ordered]) if ordered else np.zeros((0, 3))
    return seq, coords


# ---------------------------------------------------------------------
# Fixed-width PDB parsing (for the reference pdbs transferred from
# 03_protein_mpnn/inputs_prepared/)
# ---------------------------------------------------------------------
def extract_chain_ca_pdb(pdb_text, chain_id):
    by_resnum = {}
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        this_chain = line[21].strip()
        if this_chain != chain_id:
            continue
        res_name = line[17:20].strip()
        res_seq = int(line[22:26].strip())
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        by_resnum[res_seq] = (res_name, (x, y, z))
    ordered = sorted(by_resnum.items())
    seq = "".join(AA3TO1.get(resname, "X") for _, (resname, _) in ordered)
    coords = np.array([xyz for _, (_, xyz) in ordered]) if ordered else np.zeros((0, 3))
    return seq, coords


def available_pdb_chains(pdb_text):
    chains = set()
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            chains.add(line[21].strip())
    return sorted(chains)


# ---------------------------------------------------------------------
# Matchmaker-style RMSD (unchanged from the reference project)
# ---------------------------------------------------------------------
def sequence_align_correspondence(seq_a, seq_b):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(seq_a, seq_b)[0]
    pairs = []
    for (a_start, a_end), (b_start, b_end) in zip(*alignment.aligned):
        for offset in range(a_end - a_start):
            pairs.append((a_start + offset, b_start + offset))
    return pairs


def kabsch_rmsd_full(coords_a, coords_b):
    a = coords_a - coords_a.mean(axis=0)
    b = coords_b - coords_b.mean(axis=0)
    H = a.T @ b
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    a_rot = a @ R.T
    diff = a_rot - b
    return float(np.sqrt((diff ** 2).sum(axis=1).mean())), a_rot, b


def matchmaker_rmsd_multichain(chain_pairs, prune_cutoff=2.0, max_iterations=5):
    all_coords_pred, all_coords_ref = [], []
    for seq_pred, coords_pred, seq_ref, coords_ref in chain_pairs:
        if len(seq_pred) == 0 or len(seq_ref) == 0:
            continue
        pairs = sequence_align_correspondence(seq_pred, seq_ref)
        for p_idx, r_idx in pairs:
            all_coords_pred.append(coords_pred[p_idx])
            all_coords_ref.append(coords_ref[r_idx])
    if len(all_coords_pred) < 3:
        return None, 0, len(all_coords_pred)
    coords_pred_arr = np.array(all_coords_pred)
    coords_ref_arr = np.array(all_coords_ref)
    keep = np.ones(len(coords_pred_arr), dtype=bool)
    rmsd = None
    for _ in range(max_iterations):
        cp, cr = coords_pred_arr[keep], coords_ref_arr[keep]
        rmsd, cp_rot, cr_centered = kabsch_rmsd_full(cp, cr)
        per_pair_dist = np.linalg.norm(cp_rot - cr_centered, axis=1)
        new_keep_local = per_pair_dist <= prune_cutoff
        if new_keep_local.all() or new_keep_local.sum() < 3:
            break
        keep_indices = np.where(keep)[0][new_keep_local]
        keep = np.zeros(len(coords_pred_arr), dtype=bool)
        keep[keep_indices] = True
    return rmsd, int(keep.sum()), len(coords_pred_arr)


# ---------------------------------------------------------------------
# AF3 output discovery
# ---------------------------------------------------------------------
def find_model_cifs(outputs_dir, job_name):
    """<outputs_dir>/<job_name>/seed-*_sample-*/<job_name>_seed-*_sample-*_model.cif
    -- confirmed real AF3 3.0.1 naming (every file inside is prefixed with
    the job name, NOT a bare "model.cif" the way the reference project's
    output looked). Returns [(sample_index, cif_path), ...] sorted by
    sample index. Falls back to any "*_model.cif" in the sample dir (with
    a warning) if the exact prefixed name isn't found, in case the naming
    drifts again on a future AF3 version."""
    job_dir = outputs_dir / job_name
    results = []
    for sample_dir in sorted(job_dir.glob("seed-*_sample-*")):
        expected = sample_dir / f"{job_name}_{sample_dir.name}_model.cif"
        if expected.exists():
            cif_path = expected
        else:
            candidates = sorted(sample_dir.glob("*_model.cif"))
            if not candidates:
                continue
            if len(candidates) > 1:
                print(f"    [WARNING] {sample_dir}: multiple *_model.cif found, using {candidates[0].name}: "
                      f"{[c.name for c in candidates]}")
            print(f"    [NOTE] {sample_dir}: expected {expected.name} not found, using {candidates[0].name} instead")
            cif_path = candidates[0]
        sample_idx = int(sample_dir.name.rsplit("sample-", 1)[1])
        results.append((sample_idx, cif_path))
    results.sort(key=lambda x: x[0])
    return results


def load_ranking_scores(outputs_dir, job_name):
    """<outputs_dir>/<job_name>/<job_name>_ranking_scores.csv -- confirmed
    real AF3 3.0.1 naming (prefixed with the job name, not a bare
    "ranking_scores.csv"). Returns {sample_index: ranking_score}, or {} if
    the file's missing/unreadable/malformed -- callers fall back to
    sample-0 with a warning in that case, never crash."""
    job_dir = outputs_dir / job_name
    ranking_path = job_dir / f"{job_name}_ranking_scores.csv"
    if not ranking_path.exists():
        candidates = sorted(job_dir.glob("*ranking_scores.csv"))
        if not candidates:
            return {}
        print(f"    [NOTE] {job_dir}: expected {ranking_path.name} not found, using {candidates[0].name} instead")
        ranking_path = candidates[0]
    scores = {}
    with open(ranking_path) as f:
        reader = csv.DictReader(f)
        print(f"    [NOTE] {ranking_path.name} columns: {reader.fieldnames}")
        for row in reader:
            try:
                scores[int(row["sample"])] = float(row["ranking_score"])
            except (KeyError, ValueError, TypeError) as exc:
                print(f"    [WARNING] {ranking_path.name}: couldn't parse row {row}: {exc}")
    return scores


def pick_top_ranked_sample(sample_cifs, ranking_scores):
    """Higher ranking_score = better (AF3 convention). Falls back to the
    lowest sample index, with a warning, if ranking_scores.csv is
    missing/incomplete."""
    if not sample_cifs:
        return None
    scored = [(idx, path) for idx, path in sample_cifs if idx in ranking_scores]
    if not scored:
        print(f"    [WARNING] no ranking_scores.csv entries matched the found samples -- "
              f"falling back to sample-{sample_cifs[0][0]}")
        return sample_cifs[0]
    return max(scored, key=lambda entry: ranking_scores[entry[0]])


# ---------------------------------------------------------------------
# Results table -- same load-whole/write-whole CSV pattern as every
# other tracking table in this project, keyed on sequence_id.
# ---------------------------------------------------------------------
def load_table(path):
    by_sequence_id = {}
    if path.exists() and path.stat().st_size > 0:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                by_sequence_id[row["sequence_id"]] = row
    return by_sequence_id


def save_table(by_sequence_id, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_sorted = sorted(by_sequence_id.values(), key=lambda row: (row["protein_id"], row["sequence_id"]))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        for row in rows_sorted:
            writer.writerow({field: row.get(field, "") for field in RESULT_FIELDS})


# ---------------------------------------------------------------------
# outputs_clean archive merge -- same pattern as every other cumulative
# archive writer in this project.
# ---------------------------------------------------------------------
def read_archive_members(path):
    members = {}
    if path.exists():
        with tarfile.open(path, "r:gz") as archive:
            for info in archive.getmembers():
                if info.isfile():
                    handle = archive.extractfile(info)
                    members[info.name] = handle.read() if handle else b""
    return members


def write_merged_archive(archive_path, existing, new_members):
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


# ---------------------------------------------------------------------
# Discovery of sequences from this stage's own inputs/ tars -- no
# manifest.json needed, (experiment, group, protein_id, sequence_id) are
# all already known from which tar/member each sequence came from.
# ---------------------------------------------------------------------
def discover_sequences(inputs_dir):
    sequences = []
    if not inputs_dir.exists():
        return sequences
    for experiment_dir in sorted(inputs_dir.iterdir()):
        if not experiment_dir.is_dir() or experiment_dir.name == "reference":
            continue
        experiment = experiment_dir.name
        for tar_path in sorted(experiment_dir.glob("*.tar.gz")):
            group = tar_path.name[: -len(".tar.gz")]
            with tarfile.open(tar_path, "r:gz") as archive:
                for name in archive.getnames():
                    if "/" not in name or not name.endswith(".fa"):
                        continue
                    protein_id, fname = name.split("/", 1)
                    sequence_id = fname[: -len(".fa")]
                    sequences.append({
                        "protein_id": protein_id,
                        "experiment": experiment,
                        "group": group,
                        "sequence_id": sequence_id,
                    })
    return sequences


def evaluate_job(job_name, chain_pairs_fn):
    """Returns (rmsd, matched, total, ranking_score) or None if the job's
    output can't be evaluated (no samples / unparsable cif)."""
    sample_cifs = find_model_cifs(OUTPUTS_DIR, job_name)
    if not sample_cifs:
        print(f"    [WARNING] {job_name}: no AF3 samples found under {OUTPUTS_DIR / job_name}")
        return None
    ranking_scores = load_ranking_scores(OUTPUTS_DIR, job_name)
    picked = pick_top_ranked_sample(sample_cifs, ranking_scores)
    if picked is None:
        return None
    sample_idx, cif_path = picked
    try:
        rows = parse_mmcif_atom_site(cif_path.read_text())
    except ValueError as exc:
        print(f"    [WARNING] {job_name} sample-{sample_idx}: {exc}")
        return None
    rmsd, matched, total = chain_pairs_fn(rows)
    return rmsd, matched, total, ranking_scores.get(sample_idx)


def run_compute_rmsd():
    table = load_table(RESULTS_TABLE_PATH)
    already_done = set(table.keys())

    sequences = discover_sequences(INPUTS_DIR)
    print(f"[1] {len(sequences)} sequence(s) found across this stage's inputs/ tars")
    print(f"    already recorded (in {RESULTS_TABLE_PATH.name}): {len(already_done)}")

    reference_cache = {}
    new_rows = {}
    archives_to_write = {}  # (destination, experiment, group) -> {member_name: bytes}

    total_pending_af3 = 0
    total_missing_reference = 0
    total_evaluated = 0
    total_passed = 0

    for entry in sequences:
        sequence_id = entry["sequence_id"]
        if sequence_id in already_done:
            continue
        protein_id = entry["protein_id"]
        experiment = entry["experiment"]
        group = entry["group"]

        monomer_job = f"{sequence_id}_monomer"
        dimer_job = f"{sequence_id}_dimer"
        monomer_dir = OUTPUTS_DIR / monomer_job
        dimer_dir = OUTPUTS_DIR / dimer_job
        if not monomer_dir.exists() or not dimer_dir.exists():
            total_pending_af3 += 1
            continue

        if protein_id not in reference_cache:
            ref_path = REFERENCE_DIR / f"{protein_id}.pdb"
            if not ref_path.exists():
                print(f"  [WARNING] {protein_id}: no reference pdb at {ref_path} -- "
                      f"run 15082026_03_transfer_reference.py -- skipping {sequence_id}")
                total_missing_reference += 1
                continue
            ref_text = ref_path.read_text()
            ref_a = extract_chain_ca_pdb(ref_text, "A")
            ref_b = extract_chain_ca_pdb(ref_text, "B")
            if not ref_a[0] or not ref_b[0]:
                print(f"  [WARNING] {protein_id}: reference pdb has no chain A and/or B CA atoms -- "
                      f"chains actually present: {available_pdb_chains(ref_text)}")
            reference_cache[protein_id] = {"A": ref_a, "B": ref_b}
        ref = reference_cache[protein_id]

        print(f"\n--- {sequence_id} ({protein_id}, {experiment}/{group}) ---")

        def monomer_chain_pairs(rows, ref=ref):
            pred_seq, pred_coords = extract_chain_ca_mmcif(rows, "A")
            return matchmaker_rmsd_multichain([(pred_seq, pred_coords, ref["A"][0], ref["A"][1])])

        def dimer_chain_pairs(rows, ref=ref):
            pred_seq_a, pred_coords_a = extract_chain_ca_mmcif(rows, "A")
            pred_seq_b, pred_coords_b = extract_chain_ca_mmcif(rows, "B")
            return matchmaker_rmsd_multichain([
                (pred_seq_a, pred_coords_a, ref["A"][0], ref["A"][1]),
                (pred_seq_b, pred_coords_b, ref["B"][0], ref["B"][1]),
            ])

        monomer_result = evaluate_job(monomer_job, monomer_chain_pairs)
        dimer_result = evaluate_job(dimer_job, dimer_chain_pairs)
        if monomer_result is None or dimer_result is None:
            print(f"    [WARNING] {sequence_id}: could not evaluate monomer and/or dimer -- leaving for a later run")
            continue

        m_rmsd, m_matched, m_total, m_rank = monomer_result
        d_rmsd, d_matched, d_total, d_rank = dimer_result
        m_frac = (m_matched / m_total) if m_total else 0.0
        d_frac = (d_matched / d_total) if d_total else 0.0
        monomer_pass = (m_rmsd is not None) and (m_rmsd < RMSD_THRESHOLD_MONOMER) and (m_frac > MATCHED_FRACTION_THRESHOLD_MONOMER)
        dimer_pass = (d_rmsd is not None) and (d_rmsd < RMSD_THRESHOLD_DIMER) and (d_frac > MATCHED_FRACTION_THRESHOLD_DIMER)
        passed = monomer_pass and dimer_pass
        destination = "passed" if passed else "rejected"

        print(f"    monomer: RMSD={m_rmsd}, matched {m_matched}/{m_total} ({m_frac:.2%}), pass={monomer_pass}")
        print(f"    dimer:   RMSD={d_rmsd}, matched {d_matched}/{d_total} ({d_frac:.2%}), pass={dimer_pass}")
        print(f"    -> {destination}")

        # queue this sequence's full AF3 output (both jobs) for archiving
        key = (destination, experiment, group)
        bucket = archives_to_write.setdefault(key, {})
        for job_name in (monomer_job, dimer_job):
            job_dir = OUTPUTS_DIR / job_name
            for file_path in job_dir.rglob("*"):
                if file_path.is_file():
                    rel = file_path.relative_to(OUTPUTS_DIR)
                    bucket[f"{protein_id}/{rel}"] = file_path.read_bytes()

        new_rows[sequence_id] = {
            "protein_id": protein_id,
            "experiment": experiment,
            "group": group,
            "sequence_id": sequence_id,
            "monomer_rmsd": f"{m_rmsd:.4f}" if m_rmsd is not None else "",
            "monomer_matched_pairs": m_matched,
            "monomer_total_aligned_pairs": m_total,
            "monomer_matched_fraction": f"{m_frac:.4f}",
            "monomer_ranking_score": f"{m_rank:.4f}" if m_rank is not None else "",
            "monomer_pass": monomer_pass,
            "dimer_rmsd": f"{d_rmsd:.4f}" if d_rmsd is not None else "",
            "dimer_matched_pairs": d_matched,
            "dimer_total_aligned_pairs": d_total,
            "dimer_matched_fraction": f"{d_frac:.4f}",
            "dimer_ranking_score": f"{d_rank:.4f}" if d_rank is not None else "",
            "dimer_pass": dimer_pass,
            "passed": passed,
            "destination": destination,
        }
        total_evaluated += 1
        if passed:
            total_passed += 1

    # write archives
    for (destination, experiment, group), new_members in archives_to_write.items():
        archive_path = OUTPUTS_CLEAN_DIR / destination / experiment / f"{group}.tar.gz"
        existing = read_archive_members(archive_path)
        write_merged_archive(archive_path, existing, new_members)
        print(f"  archived {len(new_members)} file(s) -> {archive_path}")

    if new_rows:
        table.update(new_rows)
        save_table(table, RESULTS_TABLE_PATH)

    print(f"\nRun complete. Newly evaluated: {total_evaluated} ({total_passed} passed). "
          f"Pending AF3 completion: {total_pending_af3}. Missing reference: {total_missing_reference}.")
    print(f"  results table: {RESULTS_TABLE_PATH} ({len(table)} row(s) total)")
    print(f"  outputs_clean: {OUTPUTS_CLEAN_DIR}")


if __name__ == "__main__":
    run_compute_rmsd()
