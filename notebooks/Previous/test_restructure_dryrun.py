import sys
import numpy as np

# ── Mock BioPython structure objects ────────────────────────────────────────

class MockAtom:
    def __init__(self, coord):
        self._coord = np.array(coord, dtype=float)
    def get_vector(self):
        return self
    def get_array(self):
        return self._coord

class MockResidue:
    def __init__(self, rn, coord):
        self._id = (' ', rn, ' ')
        self._coord = coord
    def get_id(self):
        return self._id
    def __contains__(self, key):
        return key == 'CA'
    def __getitem__(self, key):
        if key == 'CA':
            return MockAtom(self._coord)
        raise KeyError(key)

class MockChain:
    def __init__(self, chain_id, residues):
        self._chain_id = chain_id
        self._residues = {r.get_id(): r for r in residues}
        self._list = residues
    def get_id(self):
        return self._chain_id
    def __getitem__(self, key):
        return self._residues[key]
    def __iter__(self):
        return iter(self._list)

class MockModel:
    def __init__(self, chain):
        self._chain = chain
    def __iter__(self):
        return iter([self._chain])
    def __getitem__(self, key):
        if key == self._chain.get_id():
            return self._chain
        raise KeyError(key)

# ── Mock data ────────────────────────────────────────────────────────────
# Q26966 has Met at pos 1, absent from all crystal constructs, so
# PDB residue N maps to Q26966 canonical position N+1 throughout.

n_residues = 620
pdb_res_nums_mock = list(range(1, n_residues + 1))
mapping_mock = {pdb_rn: pdb_rn + 1 for pdb_rn in pdb_res_nums_mock}  # pdb_rn -> canonical

# Build uniprot_seq programmatically so that canonical positions implied by
# mapping_mock actually carry the intended residue identities:
#   PDB 119 (canonical 120) = Y
#   PDB 312 (canonical 313) = W
#   PDB 342 (canonical 343) = Y   (bonus check)
uniprot_len = n_residues + 1  # +1 for the extra Met at canonical position 1
uniprot_seq_list = ['A'] * uniprot_len
uniprot_seq_list[0] = 'M'                 # canonical pos 1
uniprot_seq_list[mapping_mock[119] - 1] = 'Y'   # canonical 120
uniprot_seq_list[mapping_mock[312] - 1] = 'W'   # canonical 313
uniprot_seq_list[mapping_mock[342] - 1] = 'Y'   # canonical 343
uniprot_seq = ''.join(uniprot_seq_list)

# PDB sequence: put matching identities at the corresponding PDB (0-indexed) slots
pdb_seq_list = ['A'] * n_residues
pdb_seq_list[118] = 'Y'  # PDB res 119 (0-indexed 118)
pdb_seq_list[311] = 'W'  # PDB res 312 (0-indexed 311)
pdb_seq_list[341] = 'Y'  # PDB res 342 (0-indexed 341)
pdb_seq = ''.join(pdb_seq_list)

pdb_seed_nums = [119, 311, 312]
SPHERE_RADIUS = 18.0

REF_STRUCTURE = '1MS3'
REF_CHAIN = 'A'

# Build the reference chain geometry:
#   seed residues (PDB 119, 311, 312) at the origin
#   50 other residues within 18A
#   20 residues outside 18A
rng = np.random.default_rng(0)
used_nums = {119, 311, 312}
residues = [
    MockResidue(119, (0.0, 0.0, 0.0)),
    MockResidue(311, (0.0, 0.0, 0.0)),
    MockResidue(312, (0.0, 0.0, 0.0)),
]

next_rn = 1
def _alloc_rn():
    global next_rn
    while next_rn in used_nums:
        next_rn += 1
    rn = next_rn
    used_nums.add(rn)
    next_rn += 1
    return rn

for _ in range(50):
    rn = _alloc_rn()
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    dist = rng.uniform(0.5, 17.9)
    residues.append(MockResidue(rn, tuple(direction * dist)))

for _ in range(20):
    rn = _alloc_rn()
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    dist = rng.uniform(18.1, 40.0)
    residues.append(MockResidue(rn, tuple(direction * dist)))

ref_chain_obj_mock = MockChain(REF_CHAIN, residues)
ref_model_mock = MockModel(ref_chain_obj_mock)

structures = {
    '1MS3': {
        'uniprot_mapping': mapping_mock,
        'res_nums': pdb_res_nums_mock,
        'pdb_seq': pdb_seq,
        'structure': [ref_model_mock],   # structures[stem]['structure'][0] -> model
    }
}

# ── Cell 6c core logic (Steps 3-5 only — no GUI popups) ─────────────────────
# Steps 1-2 (PowerShell InputBox popups) are skipped in the dry run;
# pdb_seed_nums and SPHERE_RADIUS above stand in for their results.

REF_MAPPING_FWD = structures[REF_STRUCTURE]['uniprot_mapping']  # pdb_rn -> canonical

print(f'SEED RESIDUE CONVERSION — PDB numbering -> Q26966 canonical')
print('=' * 65)
print(f'Reference structure : {REF_STRUCTURE}  chain {REF_CHAIN}')
print(f'Sphere radius       : {SPHERE_RADIUS} A')
print()

SEED_RESIDUES     = []   # Q26966 canonical positions — used by all downstream cells
seed_pdb_to_canon = {}   # for reporting
seed_canon_to_pdb = {}   # reverse map for output stage

all_seeds_found = True
for pdb_rn in pdb_seed_nums:
    canonical = REF_MAPPING_FWD.get(pdb_rn, None)
    if canonical is None:
        print(f'  WARNING: PDB residue {pdb_rn} not found in alignment mapping')
        print(f'    Available PDB range: {min(REF_MAPPING_FWD.keys())} - {max(REF_MAPPING_FWD.keys())}')
        all_seeds_found = False
        continue
    canonical_aa = uniprot_seq[canonical - 1] if (canonical - 1) < len(uniprot_seq) else '?'
    ref_info = structures[REF_STRUCTURE]
    pdb_aa   = '?'
    if pdb_rn in ref_info['res_nums']:
        idx    = ref_info['res_nums'].index(pdb_rn)
        pdb_aa = ref_info['pdb_seq'][idx]
    SEED_RESIDUES.append(canonical)
    seed_pdb_to_canon[pdb_rn]  = canonical
    seed_canon_to_pdb[canonical] = pdb_rn
    match = 'OK ' if pdb_aa == canonical_aa else f'MUT({pdb_aa}!={canonical_aa})'
    print(f'  PDB {pdb_rn:>4} {pdb_aa}  ->  Q26966 pos {canonical:>4} {canonical_aa}  [{match}]')

print()
if not all_seeds_found:
    raise ValueError(
        'One or more seed residues were not found in the alignment mapping. '
        'Check that the PDB residue numbers are correct for the reference structure.'
    )

ref_model = structures[REF_STRUCTURE]['structure'][0]
if REF_CHAIN not in [c.get_id() for c in ref_model]:
    raise ValueError(f'Chain {REF_CHAIN} not found in {REF_STRUCTURE}')
ref_chain_obj = ref_model[REF_CHAIN]

seed_coords = []
for pdb_rn in pdb_seed_nums:
    try:
        res = ref_chain_obj[(' ', pdb_rn, ' ')]
        seed_coords.append(res['CA'].get_vector().get_array())
    except KeyError:
        print(f'WARNING: PDB residue {pdb_rn} not found in {REF_STRUCTURE} chain {REF_CHAIN}')

if not seed_coords:
    raise ValueError('No seed residue coordinates found. Check PDB residue numbers and chain ID.')

seed_centre = np.mean(seed_coords, axis=0)

region_res_nums = []
for r in ref_chain_obj:
    if r.get_id()[0] != ' ':
        continue
    if 'CA' not in r:
        continue
    coord = r['CA'].get_vector().get_array()
    if np.linalg.norm(coord - seed_centre) <= SPHERE_RADIUS:
        region_res_nums.append(r.get_id()[1])

region_canonical = set()
unmapped_region  = []
for pdb_rn in region_res_nums:
    canonical = REF_MAPPING_FWD.get(pdb_rn, None)
    if canonical is not None:
        region_canonical.add(canonical)
    else:
        unmapped_region.append(pdb_rn)

print(f'REGION OF INTEREST')
print('-' * 65)
print(f'Seed centre (A)         : ({seed_centre[0]:.2f}, {seed_centre[1]:.2f}, {seed_centre[2]:.2f})')
print(f'PDB residues in sphere  : {len(region_res_nums)}')
print(f'Canonical positions     : {len(region_canonical)}')
if unmapped_region:
    print(f'WARNING: {len(unmapped_region)} PDB residues could not be mapped to Q26966: {sorted(unmapped_region)}')
print()
print(f'Seed residues (PDB)       : {pdb_seed_nums}')
print(f'Seed residues (Q26966)    : {SEED_RESIDUES}')
print(f'Sphere radius             : {SPHERE_RADIUS} A')

# ── Verification checks ─────────────────────────────────────────────────────
print()
print('=' * 70)
print('DRY RUN VERIFICATION CHECKS')
print('=' * 70)

checks_passed = True

def check(n, desc, cond, actual):
    global checks_passed
    status = 'PASS' if cond else 'FAIL'
    print(f'[{status}] Check {n}: {desc} (actual={actual})')
    checks_passed &= cond

check(1, 'PDB 119 -> canonical 120, aa=Y',
      seed_pdb_to_canon.get(119) == 120 and uniprot_seq[119] == 'Y',
      (seed_pdb_to_canon.get(119), uniprot_seq[119] if len(uniprot_seq) > 119 else None))

check(2, 'PDB 312 -> canonical 313, aa=W',
      seed_pdb_to_canon.get(312) == 313 and uniprot_seq[312] == 'W',
      (seed_pdb_to_canon.get(312), uniprot_seq[312] if len(uniprot_seq) > 312 else None))

check(3, 'PDB 342 -> canonical 343, aa=Y (bonus check)',
      mapping_mock.get(342) == 343 and uniprot_seq[342] == 'Y',
      (mapping_mock.get(342), uniprot_seq[342] if len(uniprot_seq) > 342 else None))

check(4, 'region_canonical is non-empty',
      len(region_canonical) > 0, len(region_canonical))

check(5, 'SEED_RESIDUES == [120, 312, 313]',
      SEED_RESIDUES == [120, 312, 313], SEED_RESIDUES)

print()
print(f'OVERALL: {"ALL CHECKS PASSED" if checks_passed else "SOME CHECKS FAILED"}')
sys.exit(0 if checks_passed else 1)
