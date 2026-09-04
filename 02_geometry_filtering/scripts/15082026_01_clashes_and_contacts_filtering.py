import numpy as np
import os.path
import shutil
import csv
import io
import json
from datetime import datetime
import tarfile
import re
from pathlib import Path
## layer 3
"""
In this func we operate only with heavy atoms of tyrs (tyrs from the seed + tyrs determined with the output and jsons)
Kabsch transform demands next steps:
0. p -(R)-> q
1. H = p * q
2. H -(SVD)-> U, S, Vt
3. Vt.T * U.T = R
"""
CIF_COLUMNS = [
    ('record', '_atom_site.group_PDB'),
    ('element', '_atom_site.type_symbol'),
    ('label_atom', '_atom_site.label_atom_id'),
    ('label_alt', '_atom_site.label_alt_id'),
    ('label_resname', '_atom_site.label_comp_id'),
    ('label_chain', '_atom_site.label_asym_id'),
    ('label_entity', '_atom_site.label_entity_id'),
    ('label_seq', '_atom_site.label_seq_id'),
    ('ins_code', '_atom_site.pdbx_PDB_ins_code'),
    ('resi', '_atom_site.auth_seq_id'),
    ('resname', '_atom_site.auth_comp_id'),
    ('chain', '_atom_site.auth_asym_id'),
    ('atom', '_atom_site.auth_atom_id'),
    ('serial', '_atom_site.id'),
    ('b_factor', '_atom_site.B_iso_or_equiv'),
    ('occupancy', '_atom_site.occupancy'),
    ('charge', '_atom_site.pdbx_formal_charge'),
    ('is_motif', '_atom_site.is_motif_atom_with_fixed_seq'),
    ('x', '_atom_site.Cartn_x'),
    ('y', '_atom_site.Cartn_y'),
    ('z', '_atom_site.Cartn_z'),
    ('model_num', '_atom_site.pdbx_PDB_model_num'),
]
RESULT_FIELDS = [
    'protein_id',
    'experiment_name',
    'time_stamp',
    'protein_cellulose_clash_cutoff',
    'protein_protein_clash_cutoff',
    'protein_protein_contacts_cutoff',
    'rmsd_after_kabsch',
    'protein_cellulose_clashes',
    'protein_protein_clashes',
    'clashing_residues',
    'protein_protein_contacts',
    'status',
    'rejection_reason',
]
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

# ---------------------------------------------------------------------------
# CONSTANTS -- hardcoded for now, per your own decision to defer the CLI layer.
# Assumes this script lives in 02_geometry_filtering/scripts/, same
# convention as every other stage script (STAGE_ROOT = parent of scripts/).
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
STAGE_ROOT = SCRIPT_DIR.parent          # 02_geometry_filtering/
PROJECT_ROOT = STAGE_ROOT.parent        # 1cbh_clear/

INPUTS_ROOT = STAGE_ROOT / "inputs"                 # <experiment>/<group>.tar.gz, as written by the stage01->stage02 sync script
INPUTS_RAW_ROOT = STAGE_ROOT / "inputs_raw"         # scratch dir, extracted-then-deleted per protein
OUTPUTS_RAW_ROOT = STAGE_ROOT / "outputs_raw"
OUTPUTS_CLEAN_ROOT = STAGE_ROOT / "outputs_clean"   # <experiment>/{passed,rejected}/tar_groups/<group>.tar.gz, as usual
TABLE_PATH = STAGE_ROOT / "tables" / "stage_02_results.csv"

# Seed (with cellulose) is a pipeline-wide reference, not stage-specific --
# lives under the project root, not under this stage's own directory.
SEED_DIR = PROJECT_ROOT / "00_reference_structures" / "seed"
SEED_FILE = "2xseed_AB_1xcellulose_dp20_short.pdb"

CLASH_MAX = 0
CELLULOSE_CLASH_MAX = 0
CONTACT_MIN = 7
CHECKPOINT_SIZE = 100

PROTEIN_CLASH_DISTANCE_A = 1.8
CELLULOSE_CLASH_DISTANCE_A = 2.2
CONTACT_MIN_DISTANCE_A = 4.0
CONTACT_MAX_DISTANCE_A = 8.0

SEED_CHAINS = ('B', 'C')  # first = seed chain for the "directly named" generated chain, second = mirrored
def collect_dirs(root_dir: str):
    return [os.path.join(root_dir, name) for name in os.listdir(root_dir)]
def archive_direction(inputs_root: str, specific_experiment: str = None, specific_tar_gz: str = None) -> list[str]:
    inputs = os.path.normpath(inputs_root)  # helps to normalize the path
    tars = []
    if specific_tar_gz is not None:
        tars.append(os.path.join(inputs, specific_tar_gz))
        return tars
    if specific_experiment is not None:
        experiments = [os.path.join(inputs, specific_experiment)]  # os.path.join returns str by default
    else:
        experiments = collect_dirs(inputs)
    for experiment in experiments:
        tars.extend(collect_dirs(experiment))  # otherwise it would add a list in the list
    return tars
def protein_id(file_name: str) -> str:
    return (file_name.split('.'))[0]
def alignment_seed(seed_dir: str, seed_file: str) -> str:
    seed_dir = os.path.normpath(seed_dir)
    seed_path = os.path.join(seed_dir, seed_file)
    return seed_path
def parse_cif(file_path):
    clear_cif = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0:4] in ('ATOM', 'HETA'):
                line_list = line.split()
                atom = {
                    'record': line_list[0],
                    'element': line_list[1],
                    'label_atom': line_list[2],
                    'label_alt': line_list[3],
                    'label_resname': line_list[4],
                    'label_chain': line_list[5],
                    'label_entity': line_list[6],
                    'label_seq': line_list[7],
                    'ins_code': line_list[8],
                    'resi': int(line_list[9]),
                    'resname': line_list[10],
                    'chain': line_list[11],
                    'atom': line_list[12],
                    'serial': int(line_list[13]),
                    'b_factor': float(line_list[14]),
                    'occupancy': float(line_list[15]),
                    'charge': line_list[16],
                    'is_motif': line_list[17],
                    'x': float(line_list[18]),
                    'y': float(line_list[19]),
                    'z': float(line_list[20]),
                    'model_num': line_list[21],
                    }
                clear_cif.append(atom)
    clear_cif = discard_hydrogens(clear_cif)
    return clear_cif
def parse_pdb(file_path):
    clear_pdb = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0:4] in ('ATOM', 'HETA'):
                atom = {
                    'record': line[0:6].strip(),
                    'element': line[76:78].strip(),
                    'label_atom': line[12:16].strip(),
                    'label_alt': line[16].strip(),
                    'label_resname': line[17:20].strip(),
                    'label_chain': line[21].strip(),
                    'label_entity': '',
                    'label_seq': '',
                    'ins_code': line[26].strip(),
                    'resi': int(line[22:26].strip()),
                    'resname': line[17:20].strip(),
                    'chain': line[21].strip(),
                    'atom': line[12:16].strip(),
                    'serial': int(line[6:11].strip()),
                    'b_factor': float(line[60:66].strip()),
                    'occupancy': float(line[54:60].strip()),
                    'charge': '',
                    'is_motif': '',
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                    'model_num': '1',
                }
                clear_pdb.append(atom)
    clear_pdb = discard_hydrogens(clear_pdb)
    return clear_pdb
def parse_structure_file(file_path: str) -> list[dict]:
    """
    Small dispatcher: picks parse_cif or parse_pdb based on the file's
    extension, so the orchestrator doesn't need to know or care which
    format a given tar member is in.
    """
    if file_path.endswith('.cif'):
        return parse_cif(file_path)
    elif file_path.endswith('.pdb'):
        return parse_pdb(file_path)
    else:
        raise ValueError(f"Unrecognized structure file extension: {file_path}")
def write_cif(atoms: list[dict], file_path: str, protein_id: str) -> None:
    with open(file_path, 'w') as f:
        f.write(f'data_{protein_id}\n#\n')
        f.write(f'loop_\n')
        for _, field in CIF_COLUMNS:
            f.write(f'{field}\n')
        for atom in atoms:
            rows = []
            for key, _ in CIF_COLUMNS:
                value = atom.get(key, '.')
                if value == '':
                    value = '.'
                rows.append(str(value))
            f.write(' '.join(rows) + '\n')
        f.write('#\n')
def format_atom_name_pdb(atom_name: str) -> str:
    """
    Positions an atom name into the 4-character PDB atom-name field
    (cols 13-16). Simplified vs. full wwPDB atom-nomenclature rules
    (which right-justify a single-letter element symbol specifically in
    cols 13-14) -- good enough here because the only thing that ever
    reads this back is this pipeline's own parse_pdb, which just strips
    whitespace off line[12:16].
    """
    name = atom_name.strip()
    if len(name) >= 4:
        return name[:4]
    return f" {name:<3}"
def _blank_if_null(value) -> str:
    """
    CIF encodes an empty field as the literal string '.' (write_cif does
    this too, and parse_cif never translates it back), so after a
    cif -> parse_cif -> write_pdb round trip, label_alt/ins_code often
    come through as '.' rather than '' even though they're genuinely
    blank. A literal '.' in the PDB altLoc/iCode columns is non-standard
    and confuses some viewers, so treat '.' the same as '' -> a blank
    column, same as real "no altloc / no insertion code" PDB output.
    """
    text = str(value or '').strip()
    if text in ('', '.'):
        return ' '
    return text[:1]
def _pdb_ter_line(last_atom: dict) -> str:
    resname = str(last_atom.get('resname', ''))[:3].rjust(3)
    chain_id = (last_atom.get('chain', '') or '')[:1] or ' '
    resi = int(last_atom['resi'])
    return f"TER{'':>8} {resname} {chain_id}{resi:>4}\n"
def write_pdb(atoms: list[dict], file_path: str, protein_id: str) -> None:
    """
    Writes atoms as a PDB file: fixed-column ATOM/HETATM records at the
    same column offsets parse_pdb reads (record[0:6], atom[12:16],
    altloc[16], resname[17:20], chain[21], resi[22:26], icode[26],
    x/y/z[30:38 / 38:46 / 46:54], occupancy[54:60], b_factor[60:66],
    element[76:78]), a TER record whenever the chain changes, and a
    trailing END -- round-trips through this pipeline's own parse_pdb
    and is readable by normal structure viewers.
    Same non-renumbering behavior as the old write_cif: atom['serial']
    is written as-is (mod 100000 for the PDB field width), not
    reassigned, even though merge_structures combines atoms from
    multiple original sources whose serials may overlap.
    """
    with open(file_path, 'w') as f:
        f.write(f'HEADER    {protein_id}\n')
        previous_chain = None
        previous_atom = None
        for atom in atoms:
            chain = atom.get('chain', '') or ''
            if previous_chain is not None and chain != previous_chain:
                f.write(_pdb_ter_line(previous_atom))
            record = 'HETATM' if str(atom.get('record', 'ATOM')).upper().startswith('HETA') else 'ATOM'
            name = format_atom_name_pdb(str(atom.get('atom', '')))
            altloc = _blank_if_null(atom.get('label_alt', ''))
            resname = str(atom.get('resname', ''))[:3].rjust(3)
            chain_id = chain[:1] or ' '
            resi = int(atom['resi'])
            icode = _blank_if_null(atom.get('ins_code', ''))
            x, y, z = float(atom['x']), float(atom['y']), float(atom['z'])
            occupancy = float(atom.get('occupancy') or 1.0)
            b_factor = float(atom.get('b_factor') or 0.0)
            element = str(atom.get('element', ''))[:2].rjust(2)
            serial = int(atom.get('serial') or 0) % 100000
            f.write(
                f"{record:<6}{serial:>5} {name}{altloc}{resname} {chain_id}{resi:>4}{icode}   "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{b_factor:>6.2f}          {element}\n"
            )
            previous_chain = chain
            previous_atom = atom
        if previous_atom is not None:
            f.write(_pdb_ter_line(previous_atom))
        f.write('END\n')
def discard_hydrogens(structure: list[dict]) -> list[dict]:
    cleaned_structure = []
    for row in structure:
        if row['element'] not in ('H', 'D'):
            cleaned_structure.append(row)
    return cleaned_structure
# RETIRED: seed_tyr_coords, rfd3_json_matching, pairing_tyr_coords.
# These assumed one fixed anchor set (chain B/C, resi 5/31/32) and one shared
# hardcoded chain rule ('C' if atom['chain']=='A' else 'B'), and matched on
# ALL heavy atoms including CD1/CD2/CE1/CE2 -- which caused two confirmed,
# separately-measured bugs: (1) Tyr ring-atom naming ambiguity inflating the
# fitted rotation, and (2) one shared transform being the wrong model for two
# independently-diffused chains. Replaced below by index_by_chain_resi +
# build_alignment_pairs, matching your project's other script
# (align_cellulose.py), restricted to TYR_AXIS_ATOMS and driven by each
# protein's own align_keys rather than a hardcoded anchor set.
TYR_AXIS_ATOMS = ("CA", "CB", "CG", "CZ", "OH")
# CD1/CD2 and CE1/CE2 deliberately excluded: Tyr's aromatic ring can have
# symmetry-equivalent 1/2 labels or ring flips that inflate a name-matched
# RMSD. These five atoms define the backbone-to-OH axis without that
# ambiguity. (Confirmed: excluding them dropped one chain's fit RMSD from
# 1.193 A to 0.427 A on real data.)
def index_by_chain_resi(atoms: list[dict]) -> dict:
    """
    Turns a flat atom list into {chain: {resi: {atom_name: atom_dict}}} for
    O(1) lookup, instead of the old nested for-loops scanning every atom
    against every other atom.
    """
    result = {}
    for atom in atoms:
        result.setdefault(atom['chain'], {}).setdefault(atom['resi'], {})[atom['atom']] = atom
    return result
def split_map_key(key: str) -> tuple[str, int]:
    """Parses a diffused_index_map key/value like 'A24' into ('A', 24)."""
    match = re.fullmatch(r"(.+?)(-?\d+)", key)
    if not match:
        raise ValueError(f"cannot parse residue key {key!r}")
    return match.group(1), int(match.group(2))
def load_json_match(json_path: str) -> dict:
    """
    Loads one protein's matching json and returns what alignment actually
    needs: the raw diffused_index_map (chain-letter-prefixed keys/values,
    e.g. {'A24': 'B29', ...}) and this protein's own anchor residue keys
    (from specification.select_exposed, e.g. 'A5,A31,A32').
    Anchor keys are read PER PROTEIN rather than hardcoded, since not every
    protein's motif includes the same anchors (confirmed on real data --
    some proteins only fix a subset of the three Tyr positions).
    """
    with open(json_path, 'r') as f:
        data = json.load(f)
    diffused_index_map = data['diffused_index_map']
    select_exposed = data.get('specification', {}).get('select_exposed', '')
    align_keys = [key.strip() for key in select_exposed.split(',') if key.strip()]
    return {'diffused_index_map': diffused_index_map, 'align_keys': align_keys}
def build_alignment_pairs(
    seed_index: dict,
    rfd3_index: dict,
    diffused_index_map: dict,
    align_keys: list[str],
    seed_chain: str,
    rfd3_chain: str,
    fit_atom_names: tuple = TYR_AXIS_ATOMS,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Builds matched (rfd3_coords, seed_coords) arrays for ONE generated chain
    against its assigned seed chain, using only fit_atom_names (default:
    TYR_AXIS_ATOMS) -- replaces rfd3_json_matching + pairing_tyr_coords.
    Uses dict lookups (index_by_chain_resi) instead of nested loops.
    """
    rfd3_coords = []
    seed_coords = []
    for key in align_keys:
        if key not in diffused_index_map:
            continue
        _, seed_resi = split_map_key(key)
        _, rfd3_resi = split_map_key(diffused_index_map[key])
        seed_residue = seed_index.get(seed_chain, {}).get(seed_resi)
        rfd3_residue = rfd3_index.get(rfd3_chain, {}).get(rfd3_resi)
        if seed_residue is None or rfd3_residue is None:
            continue
        for atom_name in fit_atom_names:
            if atom_name in seed_residue and atom_name in rfd3_residue:
                seed_coords.append(xyz(seed_residue[atom_name]))
                rfd3_coords.append(xyz(rfd3_residue[atom_name]))
    return np.array(rfd3_coords), np.array(seed_coords)
def assign_chains(rfd3_atoms: list[dict], diffused_index_map: dict, align_keys: list[str],
                   seed_chains: tuple = ('B', 'C')) -> dict:
    """
    Determines, per protein, which generated chain is 'directly' named by
    this protein's own align_keys in diffused_index_map, and assigns the
    other generated chain to the other seed chain -- replacing the old
    hardcoded 'C' if atom['chain']=='A' else 'B' rule (confirmed worse than
    the original on real data, and non-general besides).
    Returns {generated_chain: seed_chain}.
    """
    protein_chains = sorted(set(atom['chain'] for atom in rfd3_atoms))
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
def xyz(atom) -> np.ndarray:
    try:
        return np.array([atom['x'], atom['y'], atom['z']])
    except KeyError as exc:
        raise KeyError(f"Atom missing coordinate {exc}: {atom}") from exc
def alignment_vector(rfd3_tyr_coords, seed_tyr_coords) -> tuple[np.ndarray, np.ndarray, float]:
    assert rfd3_tyr_coords.shape == seed_tyr_coords.shape
    if rfd3_tyr_coords.shape[0] < 3:
        raise ValueError(
            f"need at least 3 fitting atoms, got {rfd3_tyr_coords.shape[0]}"
        )
    centroid_rfd3_tyr_coords = np.mean(rfd3_tyr_coords, axis=0)
    centroid_seed_tyr_coords = np.mean(seed_tyr_coords, axis=0)
    rfd3_centered = rfd3_tyr_coords - centroid_rfd3_tyr_coords
    seed_centered = seed_tyr_coords - centroid_seed_tyr_coords
    covar_matrix = np.dot(rfd3_centered.T, seed_centered)
    if np.linalg.matrix_rank(covar_matrix) < 2:
        raise ValueError(
            f"fitting points are collinear or degenerate "
            f"(rank={np.linalg.matrix_rank(covar_matrix)}); cannot determine rotation"
        )
    U, S, Vt = np.linalg.svd(covar_matrix)
    if np.linalg.det(np.dot(Vt.T, U.T)) < 0.0:
        Vt[-1, :] *= -1.0
    R = np.dot(Vt.T, U.T)
    t = centroid_seed_tyr_coords - np.dot(R, centroid_rfd3_tyr_coords)
    rmsd = np.sqrt(np.sum(np.square(np.dot(rfd3_centered, R.T) - seed_centered)) / seed_centered.shape[0])
    return R, t, rmsd
def protein_alignment(rfd3_atoms, R, t) -> list[dict]:
    rfd3_transposed_protein = []
    for atom in rfd3_atoms:
        coords = xyz(atom)
        new_coords = R @ coords + t
        new_atom = atom.copy()
        new_atom.update({'x': new_coords[0], 'y': new_coords[1], 'z': new_coords[2]})
        rfd3_transposed_protein.append(new_atom)
    return rfd3_transposed_protein
def merge_structures(seed_atoms, rfd3_transposed_protein) -> list[dict]:
    rfd3_aligned_protein = list(rfd3_transposed_protein)
    for atom in seed_atoms:
        if atom['chain'] in ('K', 'L', 'M', 'N', 'O'):
            rfd3_aligned_protein.append(atom)
    return rfd3_aligned_protein
def protein_contacts(rfd3_aligned_protein: list[dict]) -> list[tuple]:
    chains = {'A': [], 'B': []}
    for atom in rfd3_aligned_protein:
        if atom['chain'] in chains and atom['atom'] == "CA":
            chains[atom['chain']].append(atom)
    for chain_id, atoms in chains.items():
        if len(atoms) == 0:
            raise ValueError(f"Chain {chain_id} is not represented in the structure")
    coords_a = np.array([xyz(atom) for atom in chains['A']])
    coords_b = np.array([xyz(atom) for atom in chains['B']])
    diff = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=2))
    interchain_contacts = []
    for i in range(dist_matrix.shape[0]):
        for j in range(dist_matrix.shape[1]):
            interchain_contacts.append((chains['A'][i]['resi'], chains['B'][j]['resi'], dist_matrix[i, j]))
    return interchain_contacts
def cellulose_contacts(rfd3_aligned_protein: list[dict]) -> list[tuple]:
    chains = {'protein': [], 'cellulose': []}
    for atom in rfd3_aligned_protein:
        if atom['chain'] in ('A', 'B'):
            chains['protein'].append(atom)
        if atom['chain'] in ('K', 'L', 'M', 'N', 'O'):
            chains['cellulose'].append(atom)
    coords_prot = np.array([xyz(atom) for atom in chains['protein']])
    coords_cell = np.array([xyz(atom) for atom in chains['cellulose']])
    diff = coords_prot[:, np.newaxis, :] - coords_cell[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=2))
    cellulose_contacts_out = []
    for i in range(dist_matrix.shape[0]):
        for j in range(dist_matrix.shape[1]):
            cellulose_contacts_out.append(
                (chains['protein'][i]['resi'], chains['protein'][i]['chain'],
                 chains['cellulose'][j]['resi'], chains['cellulose'][j]['chain'],
                 dist_matrix[i, j])
            )
    return cellulose_contacts_out
def protein_contacts_count(interchain_contacts: list[tuple]) -> int:
    num_pairs = sum(1 for pair in interchain_contacts if 4.0 <= pair[2] <= 8.0)
    return num_pairs
def protein_clashes_count(interchain_contacts: list[tuple]) -> tuple[list[tuple], int]:
    clashes = []
    c = 0
    for pair in interchain_contacts:
        if pair[2] <= 1.8:
            clashes.append((pair[0], pair[1]))
            c += 1
    return clashes, c
def cellulose_clashes_count(cellulose_contacts: list[tuple]):
    num_clashes_count = sum(1 for pair in cellulose_contacts if pair[4] <= 2.2)
    return num_clashes_count
def build_row(protein_id, experiment_name, time_stamp, rmsd, cellulose_clash_count, protein_clash_count,
              clashing_residues, contact_count, status, reasons,
              cellulose_cutoff, clash_cutoff, contact_cutoff) -> dict:
    row = {
        'protein_id': protein_id,
        'experiment_name': experiment_name,
        'time_stamp': time_stamp,
        'protein_cellulose_clash_cutoff': cellulose_cutoff,
        'protein_protein_clash_cutoff': clash_cutoff,
        'protein_protein_contacts_cutoff': contact_cutoff,
        'rmsd_after_kabsch': rmsd,
        'protein_cellulose_clashes': cellulose_clash_count,
        'protein_protein_clashes': protein_clash_count,
        'clashing_residues': str(clashing_residues),
        'protein_protein_contacts': contact_count,
        'status': status,
        'rejection_reason': str(reasons),
    }
    return row
def load_table(table_path: str) -> dict:
    rows = {}
    if os.path.exists(table_path) and os.path.getsize(table_path) > 0:
        with open(table_path, 'r') as f:
            data = csv.DictReader(f)
            for row in data:
                key = (
                    row['protein_id'],
                    float(row['protein_cellulose_clash_cutoff']),
                    float(row['protein_protein_clash_cutoff']),
                    float(row['protein_protein_contacts_cutoff'])
                )
                rows[key] = row
    return rows
def insert_row(by_id: dict, new_row: dict) -> None:
    key = (
        new_row['protein_id'],
        float(new_row['protein_cellulose_clash_cutoff']),
        float(new_row['protein_protein_clash_cutoff']),
        float(new_row['protein_protein_contacts_cutoff'])
    )
    by_id[key] = new_row
def sorted_rows(by_id: dict) -> list[dict]:
    rows = sorted(
        by_id.values(),
        key=lambda row: (
            int(row['protein_cellulose_clashes']) + int(row['protein_protein_clashes'])
            - int(row['protein_protein_contacts'])
        )
    )
    return rows
def save_table(rows: list[dict], table_path: str) -> None:
    os.makedirs(os.path.dirname(table_path), exist_ok=True)
    with open(table_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
def print_counts(by_id: dict) -> None:
    """
    PASSED/REJECTED breakdown across the whole table (all runs, not just
    this one) -- same all-time scope as "Total proteins in database"
    right next to it, and same idea as stage_01_filter.py's print_counts.
    """
    counts: dict[str, int] = {}
    for row in by_id.values():
        status = row.get('status', '')
        if status:
            counts[status] = counts.get(status, 0) + 1
    if counts:
        print("  " + ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))
def decide_status(clash_count: int, cellulose_clash_count: int, contact_count: int) -> tuple[str, list[str]]:
    reasons = []
    if clash_count > CLASH_MAX:
        reasons.append(f'more_than_{CLASH_MAX}_protein_protein_clashes')
    if cellulose_clash_count > CELLULOSE_CLASH_MAX:
        reasons.append('protein_cellulose_clashes')
    if contact_count < CONTACT_MIN:
        reasons.append(f'less_than_{CONTACT_MIN}_contacts')
    status = 'REJECTED' if reasons else 'PASSED'
    return status, reasons
def group_for_protein(protein_id: str) -> str:
    prefix = protein_id[:3]
    if prefix not in GROUP_PREFIXES:
        raise ValueError(
            f"Unknown group prefix '{prefix}' for protein_id '{protein_id}'. "
            f"Expected one of: {sorted(GROUP_PREFIXES)}"
        )
    return GROUP_PREFIXES[prefix]
def save_raw_structure(atoms: list[dict], protein_id: str, experiment_raw_dir: str,
                        outcome: str, json_source_path: str) -> tuple[str, str]:
    """
    Writes the aligned structure as PDB (write_pdb, not write_cif) and
    copies the protein's own metadata json alongside it into
    experiment_raw_dir/<outcome>/, so regroup_and_archive packs
    protein_id.json + protein_id.pdb together into each group tar.gz --
    matching the json+structure pairing every other stage in this
    pipeline (stage 01's filter, and the stage01->stage02 sync script)
    already relies on.
    """
    target_dir = os.path.join(experiment_raw_dir, outcome)
    os.makedirs(target_dir, exist_ok=True)
    pdb_path = os.path.join(target_dir, f"{protein_id}.pdb")
    write_pdb(atoms, pdb_path, protein_id)
    json_dest_path = os.path.join(target_dir, f"{protein_id}.json")
    shutil.copyfile(json_source_path, json_dest_path)
    return pdb_path, json_dest_path


def read_archive_members(archive_path: str) -> list:
    members = []
    with tarfile.open(archive_path, 'r:gz') as archive:
        for member_info in archive.getmembers():
            if not member_info.isfile():
                continue
            extracted = archive.extractfile(member_info)
            if extracted is None:
                raise RuntimeError(f"could not read member {member_info.name!r} from {archive_path}")
            members.append((member_info, extracted.read()))
    return members
def regroup_and_archive(experiment_raw_dir: str, experiment_clean_dir: str) -> None:
    """
    Cumulative across runs: MERGES this run's newly-routed loose files
    into each <group>.tar.gz rather than overwriting it -- reads any
    members already in the archive, adds this run's files, writes the
    combined result to a temp file, verifies it's complete, and
    atomically swaps it in (os.replace) before deleting the loose
    originals. Same decompress-append-recompress approach the sync
    scripts (15082026_02_transfer_to_03.py) already use on the
    destination side, applied here at the source so outputs_clean/
    itself is a complete history, not just this run's batch.
    """
    for outcome in ('passed', 'rejected'):
        raw_outcome_dir = os.path.join(experiment_raw_dir, outcome)
        if not os.path.isdir(raw_outcome_dir):
            continue
        loose_files = os.listdir(raw_outcome_dir)
        if not loose_files:
            continue
        clean_outcome_dir = os.path.join(experiment_clean_dir, outcome)
        os.makedirs(clean_outcome_dir, exist_ok=True)
        files_by_group = {}
        for file_name in loose_files:
            pid = protein_id(file_name)
            group = group_for_protein(pid)
            files_by_group.setdefault(group, []).append(file_name)
        for group, file_names in files_by_group.items():
            archive_path = os.path.join(clean_outcome_dir, f"{group}.tar.gz")

            existing_members = []
            if os.path.exists(archive_path):
                existing_members = read_archive_members(archive_path)
            existing_names = {member_info.name for member_info, _ in existing_members}

            new_names = set(file_names)
            collision = existing_names & new_names
            if collision:
                raise RuntimeError(
                    f"{archive_path}: these members are already in the archive but were "
                    f"about to be re-added: {sorted(collision)}"
                )

            temp_path = os.path.join(clean_outcome_dir, f".{group}.tar.gz.tmp{os.getpid()}")
            with tarfile.open(temp_path, 'w:gz') as archive:
                for member_info, payload in existing_members:
                    archive.addfile(member_info, io.BytesIO(payload))
                for file_name in file_names:
                    raw_path = os.path.join(raw_outcome_dir, file_name)
                    archive.add(raw_path, arcname=file_name)

            with tarfile.open(temp_path, 'r:gz') as archive:
                archived_names = set(archive.getnames())
            expected_names = existing_names | new_names
            missing = expected_names - archived_names
            if missing:
                os.remove(temp_path)
                raise RuntimeError(f"Verification failed for {archive_path}: missing {sorted(missing)}")

            os.replace(temp_path, archive_path)

            for file_name in file_names:
                os.remove(os.path.join(raw_outcome_dir, file_name))
# ---------------------------------------------------------------------------
# New helpers needed to wire everything together
# ---------------------------------------------------------------------------
def _fmt_cutoff(value: float) -> str:
    """Formats a distance like 1.8 as "1p8" so folder names avoid literal dots."""
    return f"{value:.1f}".replace(".", "p")

def output_experiment_name(input_name: str, protein_clash_distance: float, cellulose_clash_distance: float,
                            contact_min_distance: float, contact_max_distance: float) -> str:
    """
    Builds the cutoff-and-date-encoded output folder name, per your naming
    scheme: <input_name>_<date>_pp_<protein-protein clash distance, A>_pc_
    <protein-cellulose clash distance, A>_cont_<contact range min>_<max>
    e.g. 15082026_16k_production_4_seeds_4_sizes_20260817_pp_1p8_pc_2p2_cont_4p0_8p0
    These are the GEOMETRIC distance cutoffs that define a clash/contact
    (Angstroms) -- not the pass/reject count thresholds (CLASH_MAX,
    CELLULOSE_CLASH_MAX, CONTACT_MIN), which can change per run without
    changing what "a clash" or "a contact" physically means.
    Computed once per raw experiment folder name, reused for every tar
    inside that experiment during a given run.
    """
    date_str = datetime.now().strftime('%Y%m%d')
    return (
        f"{input_name}_{date_str}_pp_{_fmt_cutoff(protein_clash_distance)}_pc_{_fmt_cutoff(cellulose_clash_distance)}"
        f"_cont_{_fmt_cutoff(contact_min_distance)}_{_fmt_cutoff(contact_max_distance)}"
    )


def extract_member(archive: tarfile.TarFile, member_name: str, extract_dir: str) -> str:
    """
    Extracts one member from an already-open tar archive into extract_dir.
    Returns the resulting path on disk.
    """
    os.makedirs(extract_dir, exist_ok=True)
    archive.extract(member_name, path=extract_dir)
    return os.path.join(extract_dir, member_name)
# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
def run_filter(specific_experiment: str = None, specific_tar_gz: str = None, specific_file: str = None) -> None:
    by_id = load_table(TABLE_PATH)
    seed_path = alignment_seed(SEED_DIR, SEED_FILE)
    reference_atoms = parse_structure_file(seed_path)
    seed_index = index_by_chain_resi(reference_atoms)
    tar_paths = archive_direction(INPUTS_ROOT, specific_experiment, specific_tar_gz)
    print(f"Filter running with parameters: "
          f"protein-protein clash max={CLASH_MAX}, "
          f"cellulose clash max={CELLULOSE_CLASH_MAX}, "
          f"contact min={CONTACT_MIN}")
    output_names_cache = {}   # raw_experiment_name -> output_experiment_name
    touched_experiments = {}  # output_experiment_name -> experiment_raw_dir (only ones with new work)
    processed_since_checkpoint = 0
    first_id_this_batch = None
    last_id_this_batch = None
    for tar_path in tar_paths:
        raw_experiment_name = os.path.basename(os.path.dirname(tar_path))
        if raw_experiment_name not in output_names_cache:
                output_names_cache[raw_experiment_name] = output_experiment_name(
                raw_experiment_name, PROTEIN_CLASH_DISTANCE_A, CELLULOSE_CLASH_DISTANCE_A,
                CONTACT_MIN_DISTANCE_A, CONTACT_MAX_DISTANCE_A
            )
        experiment_name = output_names_cache[raw_experiment_name]
        with tarfile.open(tar_path) as archive:
            member_names = archive.getnames()
            for member_name in member_names:
                if not (member_name.endswith('.cif') or member_name.endswith('.pdb')):
                    continue  # skip .json members here; they're pulled in per-structure below
                pid = protein_id(member_name)
                key = (pid, float(CELLULOSE_CLASH_MAX), float(CLASH_MAX), float(CONTACT_MIN))
                if key in by_id:
                    print(f"{pid} skipped")
                    continue
                # -- extraction (temp-file approach, matching your output-side design) --
                structure_path = extract_member(archive, member_name, INPUTS_RAW_ROOT)
                json_member_name = f"{pid}.json"
                json_path = extract_member(archive, json_member_name, INPUTS_RAW_ROOT)
                rfd3_atoms = parse_structure_file(structure_path)
                json_data = load_json_match(json_path)
                # structure is fully parsed into memory -- safe to remove now.
                # json_path stays on disk a little longer: save_raw_structure
                # below needs to copy it alongside the aligned pdb output.
                os.remove(structure_path)
                # -- alignment: two independent per-chain fits, not one shared transform --
                diffused_index_map = json_data['diffused_index_map']
                align_keys = json_data['align_keys']
                rfd3_index = index_by_chain_resi(rfd3_atoms)
                chain_assignment = assign_chains(rfd3_atoms, diffused_index_map, align_keys, SEED_CHAINS)
                aligned_pieces = []
                chain_rmsds = {}
                for rfd3_chain, seed_chain in chain_assignment.items():
                    rfd3_coords, seed_coords = build_alignment_pairs(
                        seed_index, rfd3_index, diffused_index_map, align_keys,
                        seed_chain, rfd3_chain, TYR_AXIS_ATOMS,
                    )
                    if rfd3_coords.shape[0] < 3:
                        raise ValueError(
                            f"{pid}: only {rfd3_coords.shape[0]} fit atoms for chain "
                            f"{rfd3_chain} -> seed {seed_chain}; need at least 3"
                        )
                    R, t, chain_rmsd = alignment_vector(rfd3_coords, seed_coords)
                    chain_rmsds[rfd3_chain] = chain_rmsd
                    chain_atoms = [atom for atom in rfd3_atoms if atom['chain'] == rfd3_chain]
                    aligned_pieces.extend(protein_alignment(chain_atoms, R, t))
                rmsd = max(chain_rmsds.values())  # worst of the two per-chain fits, for the table/printout
                merged = merge_structures(reference_atoms, aligned_pieces)
                # -- contacts and clashes --
                interchain = protein_contacts(merged)
                clashing_residues, protein_clash_count = protein_clashes_count(interchain)
                contact_count = protein_contacts_count(interchain)
                cellulose_pairs = cellulose_contacts(merged)
                cellulose_clash_count = cellulose_clashes_count(cellulose_pairs)
                # -- pass/fail --
                status, reasons = decide_status(protein_clash_count, cellulose_clash_count, contact_count)
                outcome_dir_name = 'passed' if status == 'PASSED' else 'rejected'
                # -- per-protein printout --
                print(f"{pid}")
                print(f"  protein-protein clashes: {protein_clash_count} "
                      f"{'REJECTED' if protein_clash_count > CLASH_MAX else 'PASSED'}")
                print(f"  protein-cellulose clashes: {cellulose_clash_count} "
                      f"{'REJECTED' if cellulose_clash_count > CELLULOSE_CLASH_MAX else 'PASSED'}")
                print(f"  protein-protein contacts: {contact_count} "
                      f"{'REJECTED' if contact_count < CONTACT_MIN else 'PASSED'}")
                print(f"  Result: {status}")
                # -- assemble row, insert into in-memory table --
                new_row = build_row(
                    protein_id=pid,
                    experiment_name=experiment_name,
                    time_stamp=datetime.now().isoformat(),
                    rmsd=rmsd,
                    cellulose_clash_count=cellulose_clash_count,
                    protein_clash_count=protein_clash_count,
                    clashing_residues=clashing_residues,
                    contact_count=contact_count,
                    status=status,
                    reasons=reasons,
                    cellulose_cutoff=CELLULOSE_CLASH_MAX,
                    clash_cutoff=CLASH_MAX,
                    contact_cutoff=CONTACT_MIN,
                )
                insert_row(by_id, new_row)
                # -- write aligned structure (pdb) + json to the raw output holding area --
                experiment_raw_dir = os.path.join(OUTPUTS_RAW_ROOT, experiment_name)
                save_raw_structure(merged, pid, experiment_raw_dir, outcome_dir_name, json_path)
                os.remove(json_path)
                touched_experiments[experiment_name] = experiment_raw_dir
                # -- checkpoint bookkeeping --
                processed_since_checkpoint += 1
                if first_id_this_batch is None:
                    first_id_this_batch = pid
                last_id_this_batch = pid
                if processed_since_checkpoint >= CHECKPOINT_SIZE:
                    save_table(sorted_rows(by_id), TABLE_PATH)
                    print(f"Checkpoint: {processed_since_checkpoint} proteins processed this run "
                          f"({first_id_this_batch} to {last_id_this_batch}), table saved.")
                    processed_since_checkpoint = 0
                    first_id_this_batch = None
    # -- final flush of the table --
    if processed_since_checkpoint > 0:
        save_table(sorted_rows(by_id), TABLE_PATH)
        print(f"Final checkpoint: {processed_since_checkpoint} proteins processed "
              f"({first_id_this_batch} to {last_id_this_batch}), table saved.")
    # -- regroup + archive only the experiments that actually had new work --
    for experiment_name, experiment_raw_dir in touched_experiments.items():
        experiment_clean_dir = os.path.join(OUTPUTS_CLEAN_ROOT, experiment_name)
        regroup_and_archive(experiment_raw_dir, experiment_clean_dir)
    print(f"Run complete. Total proteins in database: {len(by_id)}")
    print_counts(by_id)
if __name__ == "__main__":
    run_filter()
