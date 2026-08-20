import pandas as pd
from collections import defaultdict

# Mock canonical sequence — 400 residues, all Ala by default
# except position 342 = Y (Tyr) and position 119 = Y (Tyr)
uniprot_seq = ['A'] * 400
uniprot_seq[118] = 'Y'   # position 119 (1-indexed) = Tyr
uniprot_seq[341] = 'Y'   # position 342 (1-indexed) = Tyr
uniprot_seq = ''.join(uniprot_seq)

# Mock structures dict — three structures:
#   1MS3: wild-type, Tyr at both positions -> both should be INCLUDED
#   3PJQ: Y342H mutant -> position 342 should be EXCLUDED, position 119 included
#   FAKE: Tyr->Phe conservative mutation at 119 -> should be EXCLUDED (Phe != Tyr)

structures = {
    '1MS3': {
        'uniprot_mapping': {119: 119, 342: 342, 200: 200},
        'sasa_per_res':    {119: 45.2, 342: 38.7, 200: 22.1},
        'res_nums':        [119, 200, 342],
        'pdb_seq':         'YAY',   # simplified: index 0=res119, 1=res200, 2=res342
    },
    '3PJQ': {
        'uniprot_mapping': {119: 119, 342: 342, 200: 200},
        'sasa_per_res':    {119: 43.1, 342: 41.2, 200: 21.8},
        'res_nums':        [119, 200, 342],
        'pdb_seq':         'YAH',   # H at position 342 = Y342H mutant
    },
    'FAKE': {
        'uniprot_mapping': {119: 119, 342: 342, 200: 200},
        'sasa_per_res':    {119: 44.0, 342: 37.9, 200: 23.0},
        'res_nums':        [119, 200, 342],
        'pdb_seq':         'FAY',   # F at position 119 = conservative Tyr->Phe mutation
    },
}

# --- New Cell 8 source (exact) ---

# Cell 8: Map SASA values to canonical UniProt positions
#
# For each structure, translate pdb_res_num -> uniprot_pos -> SASA value.
# Collect all SASA values per canonical position across all structures.
#
# MUTATION FILTER (added):
# Before recording a SASA value, verify that the PDB amino acid at this
# position matches the canonical UniProt amino acid. If they differ, this
# residue is an engineered mutation — its SASA is physiologically irrelevant
# and is excluded from the average. Only that position in that structure is
# excluded; the rest of the structure contributes normally.
#
# The alignment from Cell 6 provides everything needed:
#   structures[stem]['uniprot_mapping'] : {pdb_res_num -> uniprot_pos (1-indexed)}
#   structures[stem]['pdb_seq']         : single-letter AA string for the chain
#   structures[stem]['res_nums']        : list of pdb_res_num in sequence order
#   uniprot_seq                         : canonical UniProt sequence string

# uniprot_sasa_collection[uniprot_pos] = list of (stem, sasa, pdb_res_num, pdb_aa)
uniprot_sasa_collection = defaultdict(list)

# Track exclusions for reporting
mutation_exclusions = []   # (stem, pdb_res_num, canonical_pos, pdb_aa, canonical_aa)
total_residues_processed = 0
total_residues_excluded  = 0

for stem, info in structures.items():
    mapping     = info['uniprot_mapping']   # pdb_res_num -> uniprot_pos
    sasa_map    = info['sasa_per_res']       # pdb_res_num -> raw SASA
    pdb_seq     = info['pdb_seq']
    res_nums    = info['res_nums']

    for pdb_rn, up_pos in mapping.items():
        if pdb_rn not in sasa_map:
            continue

        total_residues_processed += 1

        # Get the amino acid in this PDB structure at this position
        idx = res_nums.index(pdb_rn) if pdb_rn in res_nums else -1
        if idx < 0:
            continue
        pdb_aa = pdb_seq[idx]

        # Get the canonical amino acid at this UniProt position
        canonical_aa = uniprot_seq[up_pos - 1] if (up_pos - 1) < len(uniprot_seq) else None

        # MUTATION FILTER: skip if amino acid differs from canonical
        # This catches all mutations regardless of chemical similarity —
        # even conservative substitutions (e.g. Tyr->Phe) are excluded
        # because the SASA reflects the mutant, not the canonical residue.
        if canonical_aa is not None and pdb_aa != canonical_aa:
            mutation_exclusions.append({
                'structure':    stem,
                'pdb_res_num':  pdb_rn,
                'uniprot_pos':  up_pos,
                'pdb_aa':       pdb_aa,
                'canonical_aa': canonical_aa,
            })
            total_residues_excluded += 1
            continue

        # Amino acid matches canonical — include this SASA value
        raw_sasa = sasa_map[pdb_rn]
        uniprot_sasa_collection[up_pos].append({
            'structure':  stem,
            'pdb_res_num': pdb_rn,
            'pdb_aa':      pdb_aa,
            'raw_sasa':    raw_sasa,
        })

# ── Report exclusions ──────────────────────────────────────────────────────────
print(f'MUTATION FILTER REPORT')
print('=' * 60)
print(f'Total residues processed : {total_residues_processed}')
print(f'Mutant residues excluded : {total_residues_excluded}')
print(f'Canonical positions with data: {len(uniprot_sasa_collection)}')
print()

if mutation_exclusions:
    excl_df = pd.DataFrame(mutation_exclusions)
    # Group by structure for a clean summary
    print('Excluded mutant residues by structure:')
    for stem, grp in excl_df.groupby('structure'):
        mutations = [
            f"{row['canonical_aa']}{int(row['uniprot_pos'])}{row['pdb_aa']}"
            for _, row in grp.iterrows()
        ]
        print(f'  {stem}: {", ".join(mutations)}')
        print(f'    (canonical_aa -> pdb_aa at UniProt position)')
    print()
    # Check for known engineered mutations in the TS dataset
    # Y342H (3PJQ), covalent intermediate modifications (2AH2)
    known_mutant_structures = excl_df['structure'].unique().tolist()
    print(f'Structures with at least one excluded residue: {known_mutant_structures}')
else:
    print('No mutant residues detected. All structures match the canonical sequence.')
    print('(If you expect mutations e.g. Y342H in 3PJQ, check that Cell 6 alignment ran correctly.)')

print()

# ── Normalise ──────────────────────────────────────────────────────────────────
# Use global max SASA across all INCLUDED (non-mutant) residues for comparability
all_raw = [
    entry['raw_sasa']
    for entries in uniprot_sasa_collection.values()
    for entry in entries
]
global_max_sasa = max(all_raw) if all_raw else 1.0
print(f'Global max SASA (mutant residues excluded): {global_max_sasa:.2f} A²')
print(f'(used for normalisation in Cell 9)')

# --- Verification checks (not part of Cell 8 itself) ---
import sys, re

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

check(1, 'total_residues_excluded == 2', total_residues_excluded == 2, total_residues_excluded)
check(2, 'len(mutation_exclusions) == 2', len(mutation_exclusions) == 2, len(mutation_exclusions))
check(3, 'canonical position 342 has exactly 2 entries',
      len(uniprot_sasa_collection.get(342, [])) == 2, len(uniprot_sasa_collection.get(342, [])))
check(4, 'canonical position 119 has exactly 2 entries',
      len(uniprot_sasa_collection.get(119, [])) == 2, len(uniprot_sasa_collection.get(119, [])))
check(5, 'canonical position 200 has exactly 3 entries',
      len(uniprot_sasa_collection.get(200, [])) == 3, len(uniprot_sasa_collection.get(200, [])))

mutation_strs = {
    f"{e['canonical_aa']}{e['uniprot_pos']}{e['pdb_aa']}" for e in mutation_exclusions
}
check(6, 'exclusion report format contains Y342H and Y119F',
      {'Y342H', 'Y119F'}.issubset(mutation_strs), mutation_strs)

check(7, 'global_max_sasa computed from included residues only (> 0)',
      global_max_sasa > 0, global_max_sasa)

print()
print(f'OVERALL: {"ALL CHECKS PASSED" if checks_passed else "SOME CHECKS FAILED"}')
sys.exit(0 if checks_passed else 1)
