import json, sys

NB = "/home/sardism/BasementAntibodyWorks/notebooks/01_target_analysis.ipynb"
OUT = "/home/sardism/BasementAntibodyWorks/data/dry_run/nb01_dry_run.py"

INJECT_CONFIG = '''
# === INJECTED (dry run, replaces popup Cell 2 + Cell 3) ===
from pathlib import Path
import numpy as np
primary_pdb_path = str(Path("~/BasementAntibodyWorks/structures/1S0I.pdb").expanduser())
alternate_pdb_path = str(Path("~/BasementAntibodyWorks/structures/3PJQ.pdb").expanduser())
output_dir = Path("~/BasementAntibodyWorks/data/dry_run/nb01").expanduser()
output_dir.mkdir(parents=True, exist_ok=True)
TARGET_NAME = "Trans-sialidase"
TARGET_LABEL = "ts"
TARGET_CHAIN = "A"
SEED_RESIDUES = [311, 119, 312]
SPHERE_RADIUS = 18.0
SASA_THRESHOLD = 0.20
EPITOPE_RADIUS = 15.0
# === END INJECTED ===
'''

INJECT_POOL = '''
# === INJECTED (dry run, replaces checkbox confirmation popup) ===
candidate_pool = candidate_df[~candidate_df['engineered']]['res_num'].tolist()
candidate_pool_df = candidate_df[~candidate_df['engineered']].copy()
# === END INJECTED ===
'''

data = json.load(open(NB))
cells = data['cells']
by_id = {c.get('id'): c for c in cells if c.get('cell_type') == 'code'}

SKIP_IDS = {'c1042938', 'afefadef'}   # Cell 2 (config form), Cell 3 (file popups)
TRUNCATE_ID = 'f85f355f'              # Cell 4c: cut before checkbox widget building
TRUNCATE_MARKER = "# Interactive checkbox table"

out_parts = []
out_parts.append("# AUTO-GENERATED dry-run script for 01_target_analysis.ipynb\n")
out_parts.append("# Non-popup cells extracted; popups replaced with injected values.\n\n")

for c in cells:
    if c.get('cell_type') != 'code':
        continue
    cid = c.get('id')
    src = ''.join(c.get('source', []))

    if cid == 'c4f4f599':
        out_parts.append(src)
        out_parts.append("\n\n")
        out_parts.append(INJECT_CONFIG)
        out_parts.append("\n")
        continue

    if cid in SKIP_IDS:
        continue

    if cid == TRUNCATE_ID:
        idx = src.find(TRUNCATE_MARKER)
        if idx == -1:
            print(f"WARNING: truncate marker not found in cell {cid}", file=sys.stderr)
            truncated = src
        else:
            truncated = src[:idx]
        out_parts.append(truncated)
        out_parts.append(INJECT_POOL)
        out_parts.append("\n")
        continue

    out_parts.append(f"\n# ---- cell {cid} ----\n")
    out_parts.append(src)
    out_parts.append("\n\n")

out_parts.append("\nprint('FINAL CANDIDATE POOL:', sorted(candidate_pool))\n")

with open(OUT, "w") as f:
    f.write("".join(out_parts))

print("wrote", OUT)
