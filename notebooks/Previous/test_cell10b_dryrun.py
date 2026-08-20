# --- mock dependencies for Cell 10b dry run ---
import pandas as pd
from IPython.display import display
display = print  # replace IPython display with print for terminal run

TARGET_NAME = 'Trans-sialidase (mock dry run)'

# candidate_pool_00a: list of UniProt canonical positions (Category A residues)
# Use a representative subset matching real TS positions from the project
candidate_pool_00a = [118, 119, 120, 122, 123, 202, 203, 249, 250,
                      280, 282, 309, 310, 312, 313, 314, 317, 362, 363]

# region_data: dict keyed by uniprot_pos (int)
# Each value must have: uniprot_residue (1-letter), mean_sasa (float),
# optionally dssp (str or None)
# Use real residue identities from UniProt Q26966 where known,
# otherwise use placeholder amino acids that exercise all tiers.
region_data = {
    118:  {'uniprot_residue': 'V', 'mean_sasa': 0.28, 'dssp': 'T'},   # Val, coil
    119:  {'uniprot_residue': 'Y', 'mean_sasa': 0.22, 'dssp': 'S'},   # Tyr, coil — force-include
    120:  {'uniprot_residue': 'N', 'mean_sasa': 0.35, 'dssp': ' '},   # Asn, coil, high SASA
    122:  {'uniprot_residue': 'G', 'mean_sasa': 0.40, 'dssp': 'T'},   # Gly, coil, high SASA
    123:  {'uniprot_residue': 'S', 'mean_sasa': 0.38, 'dssp': 'S'},   # Ser, coil, high SASA
    202:  {'uniprot_residue': 'D', 'mean_sasa': 0.32, 'dssp': 'H'},   # Asp, helix
    203:  {'uniprot_residue': 'R', 'mean_sasa': 0.29, 'dssp': 'H'},   # Arg, helix
    249:  {'uniprot_residue': 'Y', 'mean_sasa': 0.33, 'dssp': ' '},   # Tyr, coil, high SASA
    250:  {'uniprot_residue': 'W', 'mean_sasa': 0.25, 'dssp': 'E'},   # Trp, strand
    280:  {'uniprot_residue': 'T', 'mean_sasa': 0.42, 'dssp': 'T'},   # Thr, coil, high SASA
    282:  {'uniprot_residue': 'N', 'mean_sasa': 0.31, 'dssp': ' '},   # Asn, coil
    309:  {'uniprot_residue': 'L', 'mean_sasa': 0.18, 'dssp': 'E'},   # Leu, strand
    310:  {'uniprot_residue': 'A', 'mean_sasa': 0.22, 'dssp': 'E'},   # Ala, strand
    312:  {'uniprot_residue': 'W', 'mean_sasa': 0.20, 'dssp': 'S'},   # Trp, coil — force-include
    313:  {'uniprot_residue': 'S', 'mean_sasa': 0.44, 'dssp': ' '},   # Ser, coil, high SASA
    314:  {'uniprot_residue': 'G', 'mean_sasa': 0.38, 'dssp': 'T'},   # Gly, coil, high SASA
    317:  {'uniprot_residue': 'P', 'mean_sasa': 0.29, 'dssp': ' '},   # Pro, coil — penalised
    362:  {'uniprot_residue': 'Y', 'mean_sasa': 0.36, 'dssp': 'T'},   # Tyr, coil, high SASA
    363:  {'uniprot_residue': 'C', 'mean_sasa': 0.15, 'dssp': 'E'},   # Cys, strand — penalised
}

# --- Cell 10b source (exact) ---

# Cell 10b: Paratope contact scoring — rank candidate pool by antibody-interaction favorability
#
# Empirical scoring table derived from:
#   Reis et al. (2022) Front. Mol. Biosci. — 1425 Ab-Ag structures (SAbDab)
#   Kodchakorn et al. (2026) Sci. Reports — DockQ structural determinants
#
# Logic:
#   Each residue type is assigned a base paratope score reflecting how frequently
#   and versatilely that amino acid appears in real antibody paratopes.
#   The score is then modulated by two structural bonuses:
#     +30% if the residue is in a coil/loop region (Reis Fig 5D, Kodchakorn Fig 8:
#           epitopes enriched in coil; high-DockQ interfaces have 58.6% coil vs 44.3%)
#     +20% if mean relative SASA > 0.30 (Reis: >70% of epitope at medium/high exposure)
#
# Output:
#   contact_score per candidate residue, written into the hotspot JSON and
#   displayed as a ranked table. The Bayesian controller (Notebook 00) uses
#   contact_score as an acquisition prior weight so that combinations anchored
#   on Tier-1 residues (Tyr, Trp, Asp, Arg) are preferred over Tier-3 ones.
#
# References:
#   Reis et al. 2022 — Figure 4 (aa composition), Figure 6D (hydrophobic cluster),
#                      Figure 10 (Tyr/Ser interactions), Table 1
#   Kodchakorn et al. 2026 — Figure 8 (hydrophobicity), secondary structure analysis

# ── Paratope base scores ────────────────────────────────────────────────────────
# Tier 1: strongly favoured in real paratopes
# Tier 2: moderately favoured
# Tier 3: neutral / context-dependent
# Penalised: depleted in paratopes across the SAbDab dataset

PARATOPE_BASE_SCORES = {
    # Tier 1 — strongly favoured (Reis et al. Fig 4 + Table 1)
    'TYR': 1.00,  # >25% of all paratope residues; 80% in hydrophobic clusters;
                  # 40% form H-bonds; capable of pi-pi, pi-cation, pi-anion — jack of all trades
    'TRP': 0.85,  # ~10% paratope; aromatic; hydrophobic cluster participant; pi-stacking
    'SER': 0.75,  # ~10% paratope; cluster boundary role; H-bond donor/acceptor; small — no steric clash
    'ASP': 0.70,  # preferred charged residue over Glu; H-bond acceptor; salt bridge donor
    'ARG': 0.68,  # preferred over Lys; pi-cation capable; salt bridge
    'PHE': 0.65,  # aromatic; hydrophobic cluster; less versatile than Trp (lower H-bond)
    # Tier 2 — moderately favoured
    'GLY': 0.55,  # frequent in paratopes; small; no side-chain steric cost
    'THR': 0.50,  # frequent; polar; small
    'ASN': 0.50,  # frequent; polar; H-bond donor/acceptor
    'ALA': 0.40,  # neutral; small; low but non-zero contribution
    # Tier 3 — neutral / context-dependent
    'VAL': 0.30,
    'ILE': 0.25,
    'LEU': 0.25,
    # Penalised — depleted in paratopes (Reis Table 1)
    'GLN': 0.15,  # low occurrence
    'LYS': 0.15,  # less preferred charged residue vs Arg
    'GLU': 0.15,  # less preferred vs Asp
    'HIS': 0.10,  # scarce in paratopes
    'MET': 0.10,  # scarce
    'PRO': 0.05,  # very scarce; rigid backbone — conformational penalty in CDR loops
    'CYS': 0.05,  # very scarce; disulfide risk; disfavoured at interface
}

# Three-letter to one-letter map for lookup
THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

# ── Coil/loop secondary structure categories (DSSP) ────────────────────────────
# DSSP codes that represent loops/coils — these get the +30% bonus
# Reis Fig 5D: epitopes enriched in coils vs helix/strand
# Kodchakorn: high-DockQ interfaces 58.6% coil vs 44.3% in low-DockQ
COIL_DSSP_CODES = {' ', 'T', 'S', '-', 'C'}  # space/T/S = coil, turn, bend in DSSP

# ── Compute contact scores ──────────────────────────────────────────────────────
print(f'PARATOPE CONTACT SCORING — {TARGET_NAME}')
print('=' * 70)
print(f'Source: Reis et al. 2022 (1425 Ab-Ag structures), Kodchakorn et al. 2026')
print()

contact_score_records = []

for up_pos in sorted(candidate_pool_00a):
    if up_pos not in region_data:
        continue

    rd = region_data[up_pos]
    res_name_1 = rd['uniprot_residue']        # single-letter from UniProt
    res_name_3 = ONE_TO_THREE.get(res_name_1, 'UNK')

    # Base score from paratope frequency table
    base_score = PARATOPE_BASE_SCORES.get(res_name_3, 0.20)

    # Structural bonus 1: coil/loop secondary structure
    # Try to get DSSP code from region_data if it was stored; otherwise default to 0
    dssp_code = rd.get('dssp', None)
    is_coil   = (dssp_code in COIL_DSSP_CODES) if dssp_code is not None else None
    coil_bonus = 0.30 if is_coil else 0.0

    # Structural bonus 2: mean SASA above intermediate threshold
    # Reis: >70% of real epitope residues are at medium or high SASA (Rp=9 or Rp=100 probe)
    mean_sasa  = rd.get('mean_sasa', 0.0)
    sasa_bonus = 0.20 if mean_sasa > 0.30 else 0.0

    # Final score — capped at 1.0
    raw_score     = base_score * (1 + coil_bonus + sasa_bonus)
    contact_score = min(round(raw_score, 4), 1.0)

    # Tier label for display
    if base_score >= 0.75:
        tier = 'Tier 1 — strongly favoured'
    elif base_score >= 0.40:
        tier = 'Tier 2 — moderately favoured'
    elif base_score >= 0.20:
        tier = 'Tier 3 — neutral'
    else:
        tier = 'Penalised — depleted in paratopes'

    contact_score_records.append({
        'uniprot_pos':    up_pos,
        'residue_1':      res_name_1,
        'residue_3':      res_name_3,
        'base_score':     base_score,
        'mean_sasa':      round(mean_sasa, 4),
        'is_coil':        is_coil,
        'coil_bonus':     round(coil_bonus, 2),
        'sasa_bonus':     round(sasa_bonus, 2),
        'contact_score':  contact_score,
        'tier':           tier,
    })

contact_score_df = pd.DataFrame(contact_score_records).sort_values(
    'contact_score', ascending=False
).reset_index(drop=True)
contact_score_df['rank'] = contact_score_df.index + 1

print(f'Candidate pool size : {len(contact_score_df)} residues')
print()
print('RANKED CANDIDATE POOL BY CONTACT SCORE')
print('-' * 70)
display(contact_score_df[[
    'rank', 'uniprot_pos', 'residue_3', 'base_score',
    'mean_sasa', 'is_coil', 'contact_score', 'tier'
]].to_string(index=False))

print()
print('TIER SUMMARY')
print('-' * 40)
for tier_label in [
    'Tier 1 — strongly favoured',
    'Tier 2 — moderately favoured',
    'Tier 3 — neutral',
    'Penalised — depleted in paratopes',
]:
    subset = contact_score_df[contact_score_df['tier'] == tier_label]
    if len(subset) > 0:
        names = [f"{r['residue_3']}{int(r['uniprot_pos'])}"
                 for _, r in subset.iterrows()]
        print(f'{tier_label} ({len(subset)}): {", ".join(names)}')

# ── Warn if known high-value anchors are missing from Tier 1 ───────────────────
# For trans-sialidase: Tyr119, Tyr342 (catalytic), Trp312 are expected Tier 1
EXPECTED_TIER1 = []
for _, row in contact_score_df.iterrows():
    if row['residue_3'] in ('TYR', 'TRP') and row['contact_score'] >= 0.65:
        EXPECTED_TIER1.append(int(row['uniprot_pos']))

if EXPECTED_TIER1:
    print()
    print(f'Tier-1 aromatic anchors in pool (Tyr/Trp, score >= 0.65): {EXPECTED_TIER1}')
    print('These should be prioritised as anchor residues in the Bayesian controller.')
else:
    print()
    print('WARNING: No Tyr or Trp residues found in candidate pool with score >= 0.65.')
    print('Check that the candidate pool was correctly populated from Cell 10.')

# ── Append contact scores to the output JSON ───────────────────────────────────
# The JSON written in Cell 12 will be reloaded; we update it here if it already
# exists, or store the dict for Cell 12 to pick up.
contact_scores_dict = {
    str(int(row['uniprot_pos'])): row['contact_score']
    for _, row in contact_score_df.iterrows()
}

# Store on the module so Cell 12 can include it in the JSON output
import sys
_mod = sys.modules[__name__] if hasattr(sys.modules[__name__], '__file__') else None
# Fallback: just keep the variable in scope for Cell 12
# Cell 12 should add: 'contact_scores': contact_scores_dict to its output dict.

print()
print(f'contact_scores_dict populated: {len(contact_scores_dict)} entries')
print('Cell 12 will write this into the consensus SASA JSON for Notebook 00.')

# --- Verification checks (not part of Cell 10b itself) ---
print()
print('=' * 70)
print('DRY RUN VERIFICATION CHECKS')
print('=' * 70)

checks_passed = True

expected_cols = ['rank', 'uniprot_pos', 'residue_3', 'base_score',
                  'mean_sasa', 'is_coil', 'contact_score', 'tier']
ok = isinstance(contact_score_df, pd.DataFrame) and len(contact_score_df) > 0 \
     and all(c in contact_score_df.columns for c in expected_cols)
print(f'[{"PASS" if ok else "FAIL"}] contact_score_df non-empty with expected columns')
checks_passed &= ok

tier1_positions = set(
    int(r['uniprot_pos']) for _, r in contact_score_df.iterrows()
    if r['tier'] == 'Tier 1 — strongly favoured'
)
expected_tier1_positions = {119, 249, 362, 312}
ok = expected_tier1_positions.issubset(tier1_positions)
print(f'[{"PASS" if ok else "FAIL"}] Tyr119, Tyr249, Tyr362, Trp312 in Tier 1 '
      f'(got tier1={sorted(tier1_positions)})')
checks_passed &= ok

penalised_positions = set(
    int(r['uniprot_pos']) for _, r in contact_score_df.iterrows()
    if r['tier'] == 'Penalised — depleted in paratopes'
)
ok = {317, 363}.issubset(penalised_positions)
print(f'[{"PASS" if ok else "FAIL"}] Pro317, Cys363 in Penalised tier '
      f'(got penalised={sorted(penalised_positions)})')
checks_passed &= ok

ok = isinstance(contact_scores_dict, dict) and len(contact_scores_dict) == len(candidate_pool_00a)
print(f'[{"PASS" if ok else "FAIL"}] contact_scores_dict length == candidate_pool_00a length '
      f'({len(contact_scores_dict)} == {len(candidate_pool_00a)})')
checks_passed &= ok

ok = set(EXPECTED_TIER1) >= {119, 249, 362, 312}
print(f'[{"PASS" if ok else "FAIL"}] Tier-1 aromatic anchors warning block reports Tyr/Trp positions '
      f'(got {EXPECTED_TIER1})')
checks_passed &= ok

print()
print(f'OVERALL: {"ALL CHECKS PASSED" if checks_passed else "SOME CHECKS FAILED"}')
sys.exit(0 if checks_passed else 1)
