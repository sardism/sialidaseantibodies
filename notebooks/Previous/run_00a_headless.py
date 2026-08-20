#!/usr/bin/env python3
"""Headless, non-interactive runner for 00a_structural_analysisv3.ipynb.

The notebook has GUI popup cells (PowerShell InputBox / file dialogs via
subprocess) that block when run outside an interactive Jupyter session.
This script bypasses them by injecting the variables those popups would
normally set, then executes every other cell's source directly, in order,
against a shared namespace -- so it behaves like "Run All" minus the popups.

Run with: conda run -n BAW python run_00a_headless.py
"""

import sys
import traceback
from pathlib import Path

import nbformat

NOTEBOOK_PATH = Path('/home/sardism/BasementAntibodyWorks/notebooks/00a_structural_analysisv3.ipynb')

# ---------------------------------------------------------------------------
# Injected configuration -- replaces the values the popup/widget cells
# (Cell 2, Cell 3, Cell 5b, Cell 6c) would normally set interactively.
# ---------------------------------------------------------------------------
TARGET_NAME        = 'trans sialidase'
TARGET_LABEL       = 'ts'
UNIPROT_ACCESSION  = 'Q26966'   # kept for reference but not used in alignment
SASA_THRESHOLD     = 0.15
output_dir         = Path('/home/sardism/BasementAntibodyWorks/data/ts')
# NOTE: the task spec's pdb_dir (data/structures/ts/) does not exist on this
# machine. The actual 13 trans-sialidase crystal structures live in
# BasementAntibodyWorks/structures/ (confirmed against the prior consensus
# JSON's structure list). Using that directory instead.
pdb_dir             = Path('/home/sardism/BasementAntibodyWorks/structures')
REF_STRUCTURE       = '1MS3'
REF_CHAIN           = 'A'
SEED_RESIDUES       = [342, 119, 312]   # PDB residue numbers
SPHERE_RADIUS       = 18.0

output_dir.mkdir(parents=True, exist_ok=True)
pdb_paths = sorted(str(p) for p in pdb_dir.glob('*.pdb'))
pdb_files = pdb_paths  # alias matching the spec's naming

print('=' * 70)
print('HEADLESS CONFIG')
print('=' * 70)
print(f'TARGET_NAME       = {TARGET_NAME!r}')
print(f'TARGET_LABEL      = {TARGET_LABEL!r}')
print(f'UNIPROT_ACCESSION = {UNIPROT_ACCESSION!r}')
print(f'SASA_THRESHOLD    = {SASA_THRESHOLD}')
print(f'output_dir        = {output_dir}')
print(f'pdb_dir           = {pdb_dir}')
print(f'REF_STRUCTURE     = {REF_STRUCTURE}')
print(f'REF_CHAIN         = {REF_CHAIN}')
print(f'SEED_RESIDUES     = {SEED_RESIDUES}')
print(f'SPHERE_RADIUS     = {SPHERE_RADIUS}')
print(f'pdb_paths ({len(pdb_paths)}):')
for p in pdb_paths:
    print(f'  {p}')
print()

ns = dict(
    TARGET_NAME=TARGET_NAME,
    TARGET_LABEL=TARGET_LABEL,
    UNIPROT_ACCESSION=UNIPROT_ACCESSION,
    SASA_THRESHOLD=SASA_THRESHOLD,
    output_dir=output_dir,
    pdb_dir=pdb_dir,
    pdb_paths=pdb_paths,
    pdb_files=pdb_files,
    REF_STRUCTURE=REF_STRUCTURE,
    REF_CHAIN=REF_CHAIN,
    SEED_RESIDUES=SEED_RESIDUES,
    SPHERE_RADIUS=SPHERE_RADIUS,
)

# Force a non-interactive matplotlib backend before Cell 1 imports pyplot.
import matplotlib
matplotlib.use('Agg')
ns['__name__'] = '__main__'

# ---------------------------------------------------------------------------
# Popup/widget cells to skip entirely -- their effect is replaced by the
# injected config above. Matched by a unique substring of each cell's
# leading comment.
# ---------------------------------------------------------------------------
SKIP_MARKERS = [
    '# Cell 2: Configuration',
    '# Cell 3: Select output directory and PDB files',
    '# Cell 5b: Reference structure and region of interest selection',
]

# Cell 6c mixes a popup (seed residues + sphere radius input) with real
# computation (locating the region of interest around the seed residues).
# We keep the computation, replacing the popup-derived `pdb_seed_nums`
# with the pre-injected SEED_RESIDUES/SPHERE_RADIUS.
CELL_6C_MARKER = '# Cell 6c: Seed residue selection and region of interest definition'
CELL_6C_REPLACEMENT = '''
pdb_seed_nums = list(SEED_RESIDUES)
if not pdb_seed_nums:
    raise ValueError('No seed residues entered. Re-run this cell.')

print(f'SEED RESIDUES — PDB numbering')
print('=' * 65)
print(f'Reference structure : {REF_STRUCTURE}  chain {REF_CHAIN}')
print(f'Sphere radius       : {SPHERE_RADIUS} A')
print()

ref_info = structures[REF_STRUCTURE]
all_found = True
for pdb_rn in pdb_seed_nums:
    if pdb_rn in ref_info['res_nums']:
        idx    = ref_info['res_nums'].index(pdb_rn)
        pdb_aa = ref_info['pdb_seq'][idx]
        print(f'  PDB {pdb_rn:>4} = {pdb_aa}')
    else:
        print(f'  WARNING: PDB {pdb_rn} not found in reference structure')
        all_found = False

print()
if not all_found:
    raise ValueError('One or more seed residues not found. Check PDB residue numbers.')

ref_model     = structures[REF_STRUCTURE]['structure'][0]
ref_chain_obj = ref_model[REF_CHAIN]

seed_coords = []
for pdb_rn in pdb_seed_nums:
    try:
        res = ref_chain_obj[(' ', pdb_rn, ' ')]
        seed_coords.append(res['CA'].get_vector().get_array())
    except KeyError:
        print(f'WARNING: PDB {pdb_rn} not found in {REF_STRUCTURE} chain {REF_CHAIN}')

if not seed_coords:
    raise ValueError('No seed coordinates found.')

seed_centre = np.mean(seed_coords, axis=0)

region_res_nums = []
for r in ref_chain_obj:
    if r.get_id()[0] != ' ' or 'CA' not in r:
        continue
    coord = r['CA'].get_vector().get_array()
    if np.linalg.norm(coord - seed_centre) <= SPHERE_RADIUS:
        region_res_nums.append(r.get_id()[1])

region_canonical = set(region_res_nums)

print(f'REGION OF INTEREST')
print('-' * 65)
print(f'Seed centre (A) : ({seed_centre[0]:.2f}, {seed_centre[1]:.2f}, {seed_centre[2]:.2f})')
print(f'PDB residues in sphere: {len(region_res_nums)}')
print(f'Seed residues   : {SEED_RESIDUES}')
print(f'Sphere radius   : {SPHERE_RADIUS} A')
print()
print('All numbers are PDB residue numbers.')
print('Tyr342=342, Tyr119=119, Trp312=312 throughout.')
'''


def cell_label(src: str) -> str:
    return src.split('\n')[0][:70]


def main():
    nb = nbformat.read(NOTEBOOK_PATH, as_version=4)
    code_cells = [c for c in nb.cells if c.cell_type == 'code']

    for i, cell in enumerate(code_cells):
        src = cell.source
        label = cell_label(src)

        if any(src.startswith(marker) for marker in SKIP_MARKERS):
            print(f'--- [{i}] SKIP (popup/widget cell): {label}')
            continue

        if src.startswith(CELL_6C_MARKER):
            print(f'--- [{i}] RUN (popup stripped): {label}')
            src = CELL_6C_REPLACEMENT
        else:
            print(f'--- [{i}] RUN: {label}')

        try:
            exec(compile(src, f'<cell {i}: {label}>', 'exec'), ns)
        except Exception:
            print()
            print(f'!!! ERROR in cell [{i}]: {label}')
            traceback.print_exc()
            sys.exit(1)

    print()
    print('=' * 70)
    print('ALL CELLS EXECUTED SUCCESSFULLY')
    print('=' * 70)

    # -----------------------------------------------------------------
    # Verify the output JSON was written and has the expected shape.
    # -----------------------------------------------------------------
    out_path = output_dir / f'{TARGET_LABEL}_consensus_sasa.json'
    if not out_path.exists():
        print(f'!!! Output JSON not found: {out_path}')
        sys.exit(1)

    import json
    with open(out_path) as f:
        data = json.load(f)

    rd = data.get('region_data', {})
    checks = [
        ('candidate_pool_A has >=1 residue', len(data.get('candidate_pool_A', [])) > 0),
        ('contact_scores dict non-empty', len(data.get('contact_scores', {})) > 0),
        ('residue 342 present', '342' in rd),
        ('residue 119 present', '119' in rd),
        ('residue 312 present', '312' in rd),
        ('Tyr at 342', rd.get('342', {}).get('uniprot_residue') == 'Y'),
        ('Tyr at 119', rd.get('119', {}).get('uniprot_residue') == 'Y'),
        ('Trp at 312', rd.get('312', {}).get('uniprot_residue') == 'W'),
    ]

    print()
    print('OUTPUT JSON CHECKS')
    print('-' * 70)
    all_pass = True
    for name, ok in checks:
        print(f'  [{"PASS" if ok else "FAIL"}] {name}')
        all_pass = all_pass and ok

    print()
    if all_pass:
        print('ALL CHECKS PASSED')
        print(f'Candidate pool: {data["candidate_pool_A"]}')
        print(f'Contact scores: {data["contact_scores"]}')
    else:
        print('SOME CHECKS FAILED')
        sys.exit(1)


if __name__ == '__main__':
    main()
