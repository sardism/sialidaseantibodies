# AUTO-GENERATED dry-run script for 01_target_analysis.ipynb
# Non-popup cells extracted; popups replaced with injected values.

# Cell 1: Imports and popup helpers

import os
import sys
import json
import subprocess
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

from scipy.spatial.distance import pdist, squareform

from Bio.PDB import PDBParser, PPBuilder, PDBIO, Select, NeighborSearch

import ipywidgets as widgets
from IPython.display import display, HTML, clear_output

warnings.filterwarnings('ignore')

# --------------------------------------------------------------------------
# Popup helpers -- native Windows dialogs via PowerShell (WSL interop)
# --------------------------------------------------------------------------
def ps_browse_file(title='Select file', filter_str='PDB files|*.pdb|JSON files|*.json|All files|*.*'):
    script = f"""
Add-Type -AssemblyName System.Windows.Forms
$d = New-Object System.Windows.Forms.OpenFileDialog
$d.Title = "{title}"
$d.Filter = "{filter_str}"
$d.ShowDialog() | Out-Null
Write-Output $d.FileName
"""
    r = subprocess.run(["powershell.exe", "-NoProfile", "-Command", script], capture_output=True, text=True)
    p = r.stdout.strip()
    if not p:
        return None
    w = subprocess.run(["wslpath", p], capture_output=True, text=True)
    return w.stdout.strip() or None

def ps_browse_folder(title='Select folder'):
    script = f"""
Add-Type -AssemblyName System.Windows.Forms
$d = New-Object System.Windows.Forms.FolderBrowserDialog
$d.Description = "{title}"
$d.ShowNewFolderButton = $true
$d.ShowDialog() | Out-Null
Write-Output $d.SelectedPath
"""
    r = subprocess.run(["powershell.exe", "-NoProfile", "-Command", script], capture_output=True, text=True)
    p = r.stdout.strip()
    if not p:
        return None
    w = subprocess.run(["wslpath", p], capture_output=True, text=True)
    return w.stdout.strip() or None

def ps_multiline_input(title='Input', prompt='Enter data:', default_text='', width=700, height=550):
    safe = default_text.replace('"', '`"').replace('\n', '`r`n')
    bw, bh = width - 110, height - 80
    cw, ch = width - 210, height - 80
    script = f"""
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
$f = New-Object System.Windows.Forms.Form
$f.Text = "{title}"; $f.Size = New-Object System.Drawing.Size({width},{height})
$f.StartPosition = "CenterScreen"; $f.TopMost = $true
$lb = New-Object System.Windows.Forms.Label
$lb.Text = "{prompt}"; $lb.Location = New-Object System.Drawing.Point(10,10)
$lb.Size = New-Object System.Drawing.Size({width-30},30); $f.Controls.Add($lb)
$tb = New-Object System.Windows.Forms.TextBox
$tb.Multiline = $true; $tb.ScrollBars = "Vertical"
$tb.Location = New-Object System.Drawing.Point(10,45)
$tb.Size = New-Object System.Drawing.Size({width-30},{height-130})
$tb.Font = New-Object System.Drawing.Font("Consolas",10)
$tb.Text = "{safe}"; $f.Controls.Add($tb)
$ok = New-Object System.Windows.Forms.Button
$ok.Text = "Confirm"; $ok.Location = New-Object System.Drawing.Point({bw},{bh})
$ok.Size = New-Object System.Drawing.Size(90,32)
$ok.DialogResult = [System.Windows.Forms.DialogResult]::OK; $f.Controls.Add($ok)
$ca = New-Object System.Windows.Forms.Button
$ca.Text = "Cancel"; $ca.Location = New-Object System.Drawing.Point({cw},{ch})
$ca.Size = New-Object System.Drawing.Size(90,32)
$ca.DialogResult = [System.Windows.Forms.DialogResult]::Cancel; $f.Controls.Add($ca)
$f.AcceptButton = $ok; $f.CancelButton = $ca
if ($f.ShowDialog() -eq [System.Windows.Forms.DialogResult]::OK) {{ Write-Output $tb.Text }}
"""
    r = subprocess.run(["powershell.exe", "-NoProfile", "-Command", script], capture_output=True, text=True)
    return r.stdout.strip()

print('Imports loaded.')
print(f'BioPython PDB parser ready.')
print(f'NumPy {np.__version__}, pandas {pd.__version__}')



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


# ---- cell 5f28b09f ----
# Cell 4: Load structure and pick chain

parser = PDBParser(QUIET=True)
ppb    = PPBuilder()

structure = parser.get_structure(TARGET_LABEL, primary_pdb_path)
model     = structure[0]

print(f'{TARGET_NAME.upper()} STRUCTURE')
print('=' * 60)
print(f'File  : {Path(primary_pdb_path).name}')
print(f'Models: {len(structure)}')
print()

chain_summary = []
for chain in model:
    cid = chain.get_id()
    protein_res = [r for r in chain if r.get_id()[0] == ' ']
    hetatm_res  = [r for r in chain if r.get_id()[0] not in (' ', 'W')]
    water_res   = [r for r in chain if r.get_id()[0] == 'W']

    seq = ''
    for pp in ppb.build_peptides(chain):
        seq += str(pp.get_sequence())

    res_range = ''
    if protein_res:
        first = protein_res[0].get_id()[1]
        last  = protein_res[-1].get_id()[1]
        res_range = f'{first}–{last}'

    ligands = set(r.get_resname() for r in hetatm_res)

    print(f'Chain {cid}:')
    print(f'  Protein residues : {len(protein_res)}  ({res_range})')
    print(f'  Sequence length  : {len(seq)}')
    print(f'  HETATM residues  : {len(hetatm_res)} {ligands if ligands else ""}')
    print(f'  Waters           : {len(water_res)}')
    print()

    chain_summary.append({'id': cid, 'n_residues': len(protein_res), 'seq_len': len(seq)})

# ---- Chain picker widget ----
chain_ids   = [c['id'] for c in chain_summary]
default_idx = chain_ids.index(TARGET_CHAIN) if TARGET_CHAIN in chain_ids else 0

chain_picker = widgets.Dropdown(
    options=chain_ids,
    value=chain_ids[default_idx],
    description='Work with chain:',
    style={'description_width': 'initial'},
)

confirm_btn  = widgets.Button(description='Confirm chain', button_style='success')
chain_output = widgets.Output()

def on_confirm(_):
    global TARGET_CHAIN, working_chain
    TARGET_CHAIN  = chain_picker.value
    working_chain = model[TARGET_CHAIN]
    with chain_output:
        clear_output()
        print(f'Chain {TARGET_CHAIN} selected. {len([r for r in working_chain if r.get_id()[0]==" "])} protein residues.')

confirm_btn.on_click(on_confirm)
display(widgets.HBox([chain_picker, confirm_btn]), chain_output)

# Initialise silently so later cells work even if button is not clicked
TARGET_CHAIN  = chain_ids[default_idx]
working_chain = model[TARGET_CHAIN]



# ---- cell e85e6009 ----
# Cell 4b: Surface-exposed residues (SASA)
#
# Computes solvent-accessible surface area for every residue in the selected
# chain using BioPython's Shrake-Rupley algorithm, and records each residue's
# CA coordinates. This feeds the candidate-pool search in Cell 4c.

from Bio.PDB.SASA import ShrakeRupley
import numpy as np

SASA_THRESHOLD_VAL = SASA_THRESHOLD  # from Cell 2

sr = ShrakeRupley()
sr.compute(structure, level="R")

sasa_data = []
for residue in working_chain:
    if residue.get_id()[0] != ' ' or 'CA' not in residue:
        continue
    res_num = residue.get_id()[1]
    res_name = residue.get_resname()
    sasa = getattr(residue, 'sasa', 0.0)
    ca_coords = residue['CA'].get_vector().get_array()
    sasa_data.append({
        'res_num': res_num,
        'res_name': res_name,
        'sasa': sasa,
        'x': ca_coords[0],
        'y': ca_coords[1],
        'z': ca_coords[2],
    })

sasa_df = pd.DataFrame(sasa_data)
max_sasa = sasa_df['sasa'].max()
sasa_df['relative_sasa'] = sasa_df['sasa'] / (max_sasa + 1e-9)
sasa_df['exposed'] = sasa_df['relative_sasa'] >= SASA_THRESHOLD_VAL

print(f'SASA ANALYSIS — {TARGET_NAME} chain {TARGET_CHAIN}')
print(f'Total residues     : {len(sasa_df)}')
print(f'Surface-exposed    : {sasa_df["exposed"].sum()} ({100*sasa_df["exposed"].mean():.0f}%)')
print(f'Buried             : {(~sasa_df["exposed"]).sum()}')
print()
print('Top 20 most exposed residues:')
display(sasa_df.sort_values('sasa', ascending=False).head(20)[['res_num','res_name','sasa','relative_sasa']].round(3))


# Cell 4c: Candidate residue pool
#
# Finds all surface-exposed residues within SPHERE_RADIUS of the seed residues,
# shows an interactive checkbox table, and locks in the candidate pool on
# Confirm. This pool is what Notebook 00 (the Bayesian controller) explores.

# Parse seed residues
seed_res_nums = SEED_RESIDUES  # list of ints from Cell 2

# Get seed coordinates
seed_coords = []
for res_num in seed_res_nums:
    try:
        res = working_chain[(' ', res_num, ' ')]
        seed_coords.append(res['CA'].get_vector().get_array())
    except KeyError:
        print(f'WARNING: seed residue {res_num} not found in chain {TARGET_CHAIN}')

if not seed_coords:
    raise ValueError('No seed residues found. Check residue numbers match the PDB.')

seed_centre = np.mean(seed_coords, axis=0)

# Collect residues within sphere
in_sphere = []
for _, row in sasa_df.iterrows():
    coord = np.array([row['x'], row['y'], row['z']])
    dist = np.linalg.norm(coord - seed_centre)
    in_sphere.append(dist <= SPHERE_RADIUS)

sasa_df['in_sphere'] = in_sphere
sasa_df['dist_to_centre'] = sasa_df.apply(
    lambda r: np.linalg.norm(np.array([r['x'],r['y'],r['z']]) - seed_centre), axis=1)

# Known engineered mutation positions in 1S0I to flag
ENGINEERED_MUTATIONS = [58, 59, 495, 496, 520, 593, 597, 599]

candidate_df = sasa_df[sasa_df['in_sphere'] & sasa_df['exposed']].copy()
candidate_df = candidate_df.sort_values('sasa', ascending=False).reset_index(drop=True)
candidate_df['engineered'] = candidate_df['res_num'].isin(ENGINEERED_MUTATIONS)
candidate_df['include'] = ~candidate_df['engineered']

print(f'CANDIDATE POOL — {TARGET_NAME}')
print(f'Seed residues     : {seed_res_nums}')
print(f'Sphere radius     : {SPHERE_RADIUS} A')
print(f'SASA threshold    : {SASA_THRESHOLD_VAL}')
print(f'Candidates found  : {len(candidate_df)}')
print(f'Engineered (flagged excluded): {candidate_df["engineered"].sum()}')
print()


# === INJECTED (dry run, replaces checkbox confirmation popup) ===
candidate_pool = candidate_df[~candidate_df['engineered']]['res_num'].tolist()
candidate_pool_df = candidate_df[~candidate_df['engineered']].copy()
# === END INJECTED ===


# ---- cell f25f039d ----
# Cell 5: Epitope centre of mass
#
# The candidate pool selected in Cell 4c is what Notebook 00 (the Bayesian
# controller) explores. Here we simply fix the epitope centre of mass from the
# seed residues and expose the candidate pool for the downstream geometry and
# neighbourhood cells.

epitope_com = seed_centre.copy()
epitope_core_nums = candidate_pool
print(f'Epitope centre of mass: ({epitope_com[0]:.2f}, {epitope_com[1]:.2f}, {epitope_com[2]:.2f})')
print(f'Candidate pool size: {len(candidate_pool)} residues')
print(f'Residues: {sorted(candidate_pool)}')

# Empty core-residue table kept so the geometry cell (Cell 6) still runs;
# in this workflow geometry is computed from the seed-derived centre.
epitope_core_df = pd.DataFrame()
epitope_core_data = []



# ---- cell cc52e990 ----
# Cell 6: Epitope geometry

if not epitope_core_df.empty:
    matched = epitope_core_df if 'match' not in epitope_core_df.columns else \
              epitope_core_df[epitope_core_df['match'] == True]
else:
    matched = pd.DataFrame()

if len(matched) > 1:
    coords = matched[['x', 'y', 'z']].values

    from scipy.spatial.distance import pdist, squareform
    pairwise = squareform(pdist(coords))

    dists_from_com  = np.linalg.norm(coords - epitope_com, axis=1)
    max_spread      = dists_from_com.max()
    mean_spread     = dists_from_com.mean()

    print(f'EPITOPE GEOMETRY — {TARGET_NAME}')
    print('=' * 60)
    print(f'Centre of mass:       ({epitope_com[0]:.2f}, {epitope_com[1]:.2f}, {epitope_com[2]:.2f})')
    print(f'Max spread from COM:   {max_spread:.1f} Angstroms')
    print(f'Mean spread:           {mean_spread:.1f} Angstroms')
    print(f'Max pairwise dist:     {pairwise.max():.1f} Angstroms')
    print()

    labels = [f"{row['actual_aa']}{row['res_num']}" for _, row in matched.iterrows()]
    pw_df  = pd.DataFrame(pairwise, index=labels, columns=labels)
    print('Pairwise CA distances (Angstroms):')
    display(pw_df.round(1))
elif len(matched) == 1:
    print('Only one core residue — pairwise geometry not applicable.')
    print(f'Centre: ({epitope_com[0]:.2f}, {epitope_com[1]:.2f}, {epitope_com[2]:.2f})')
else:
    print(f'Epitope centre: ({epitope_com[0]:.2f}, {epitope_com[1]:.2f}, {epitope_com[2]:.2f})')
    print('(No core residue table — geometry computed from geometric/auto centre.)')



# ---- cell f3da996b ----
# Cell 7: Find all residues in the epitope neighbourhood
#
# Collects every residue within EPITOPE_RADIUS of the epitope centre.
# This is the full region an antibody could contact.

all_ca_atoms = [
    residue['CA']
    for residue in working_chain
    if residue.get_id()[0] == ' ' and 'CA' in residue
]

ns             = NeighborSearch(all_ca_atoms)
nearby_residues = ns.search(epitope_com.tolist(), EPITOPE_RADIUS, level='R')

neighbourhood_data = []
for residue in nearby_residues:
    if residue.get_id()[0] != ' ':
        continue
    res_num    = residue.get_id()[1]
    res_name   = residue.get_resname()
    ca_coords  = residue['CA'].get_vector().get_array()
    dist_to_com = np.linalg.norm(ca_coords - epitope_com)
    in_core    = res_num in epitope_core_nums

    neighbourhood_data.append({
        'res_num':         res_num,
        'res_name':        res_name,
        'dist_to_epitope': dist_to_com,
        'in_core':         in_core,
        'x': ca_coords[0], 'y': ca_coords[1], 'z': ca_coords[2],
    })

neighbourhood_df = pd.DataFrame(neighbourhood_data).sort_values('dist_to_epitope').reset_index(drop=True)

print(f'EPITOPE NEIGHBOURHOOD — {TARGET_NAME}')
print(f'Radius: {EPITOPE_RADIUS} A from centre')
print('=' * 60)
print(f'Total residues in neighbourhood : {len(neighbourhood_df)}')
print(f'Core epitope residues           : {neighbourhood_df["in_core"].sum()}')
print(f'Surrounding residues            : {(~neighbourhood_df["in_core"]).sum()}')
print()
display(neighbourhood_df)



# ---- cell fc9286f3 ----
# Cell 9: Compare two structures (optional)
#
# If you loaded an alternate structure in Cell 3, this cell compares
# the two structures residue by residue at the candidate-pool positions.
#
# Use cases:
#   - Wild-type vs mutant: check if the epitope is conserved
#   - Apo vs holo: see conformational changes at the binding site
#   - Active vs inactive: identify conformation-specific residues

if alternate_pdb_path:
    alt_structure = parser.get_structure(f'{TARGET_LABEL}_alt', alternate_pdb_path)
    alt_model     = alt_structure[0]

    # Use the same chain ID as the primary structure (adjust if needed)
    alt_chains = [c.get_id() for c in alt_model]
    alt_chain_id = TARGET_CHAIN if TARGET_CHAIN in alt_chains else alt_chains[0]
    alt_chain  = alt_model[alt_chain_id]

    print(f'STRUCTURE COMPARISON — {TARGET_NAME}')
    print('=' * 60)
    print(f'Primary  : {Path(primary_pdb_path).name}   chain {TARGET_CHAIN}')
    print(f'Alternate: {Path(alternate_pdb_path).name}  chain {alt_chain_id}')
    print()

    # Compare across the candidate pool residues
    compare_nums = epitope_core_nums

    print(f'Comparing {len(compare_nums)} residue position(s).')
    print()
    print(f'{"Position":<10} {"Primary":<10} {"Alternate":<12} {"Same?"}')
    print('-' * 45)

    differences = []
    for res_num in sorted(compare_nums):
        primary_aa  = '---'
        alternate_aa = '---'

        try:
            primary_aa = working_chain[(' ', res_num, ' ')].get_resname()
        except KeyError:
            pass

        try:
            alternate_aa = alt_chain[(' ', res_num, ' ')].get_resname()
        except KeyError:
            try:
                alternate_aa = alt_chain[res_num].get_resname()
            except:
                pass

        same = 'YES' if primary_aa == alternate_aa else 'DIFFERENT'
        if same == 'DIFFERENT':
            differences.append(res_num)
        print(f'{res_num:<10} {primary_aa:<10} {alternate_aa:<12} {same}')

    print()
    if differences:
        print(f'Positions that differ: {differences}')
        print()
        print('TIP: If residues differ at candidate positions, the antibody may')
        print('show differential binding between the two structural states.')
        print('If residues are conserved, the antibody should bind both states.')
    else:
        print('All compared positions are identical between the two structures.')

else:
    print('No alternate structure provided. Skipping comparison.')
    print('Select an alternate PDB in Cell 3 to enable this analysis.')



# ---- cell 7cc43c75 ----
# Cell 10: Clean and export the target PDB
#
# RFdiffusion and Boltz-2 need a clean PDB:
#   - Selected chain only
#   - Protein residues only (no waters, no ligands, no HETATM)
#   - No alternate conformations (keep A or blank)

class CleanProteinSelect(Select):
    """Select protein residues from one chain; drop waters, ligands, altlocs."""
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

    def accept_residue(self, residue):
        return residue.get_id()[0] == ' '   # space = standard amino acid

    def accept_atom(self, atom):
        altloc = atom.get_altloc()
        return altloc == ' ' or altloc == 'A'

pdb_stem     = Path(primary_pdb_path).stem
clean_pdb_path = output_dir / f'{pdb_stem}_chain{TARGET_CHAIN}_clean.pdb'

io = PDBIO()
io.set_structure(structure)
io.save(str(clean_pdb_path), CleanProteinSelect(TARGET_CHAIN))

# Verify
clean_struct    = parser.get_structure('clean', str(clean_pdb_path))
clean_chain     = clean_struct[0][TARGET_CHAIN]
clean_residues  = [r for r in clean_chain if r.get_id()[0] == ' ']

print(f'CLEANED TARGET PDB — {TARGET_NAME}')
print('=' * 60)
print(f'Output file : {clean_pdb_path}')
print(f'Chain       : {TARGET_CHAIN}')
print(f'Residues    : {len(clean_residues)}')
print(f'Range       : {clean_residues[0].get_id()[1]} — {clean_residues[-1].get_id()[1]}')
print()
print('This file is ready for RFdiffusion and Boltz-2.')



# ---- cell 5b84411b ----
# Cell 11: Save all analysis outputs

# 1. Core epitope residues CSV (only if a core table was built)
if not epitope_core_df.empty:
    core_csv = output_dir / f'{TARGET_LABEL}_epitope_core.csv'
    epitope_core_df.to_csv(core_csv, index=False)
    print(f'Core epitope CSV   : {core_csv}')

# 2. Full neighbourhood CSV
neighbourhood_csv = output_dir / f'{TARGET_LABEL}_epitope_neighbourhood.csv'
neighbourhood_df.to_csv(neighbourhood_csv, index=False)
print(f'Neighbourhood CSV  : {neighbourhood_csv}')

# 3. Candidate pool JSON (input for Notebook 00 / Notebook 02)
hotspot_json_path = output_dir / f'{TARGET_LABEL}_hotspot_residues.json'
hotspot_data = {
    'target_name': TARGET_NAME,
    'source_pdb': Path(primary_pdb_path).name,
    'chain': TARGET_CHAIN,
    'seed_residues': seed_res_nums,
    'sphere_radius': SPHERE_RADIUS,
    'sasa_threshold': SASA_THRESHOLD_VAL,
    'epitope_com': epitope_com.tolist(),
    'candidate_pool': sorted(candidate_pool),
    'n_candidates': len(candidate_pool),
    'clean_pdb': str(clean_pdb_path),
}
with open(hotspot_json_path, 'w') as f:
    json.dump(hotspot_data, f, indent=2)
print(f'Candidate pool JSON: {hotspot_json_path}')

# 4. Epitope centre of mass (plain text for tools that need raw coords)
com_file = output_dir / f'{TARGET_LABEL}_epitope_com.txt'
com_file.write_text(f'{epitope_com[0]:.3f} {epitope_com[1]:.3f} {epitope_com[2]:.3f}\n')
print(f'Epitope COM        : {com_file}')

print()
print('=' * 60)
print(f'NOTEBOOK 01 COMPLETE — {TARGET_NAME}')
print('=' * 60)
print()
print(f'  Source PDB            : {Path(primary_pdb_path).name}')
print(f'  Chain                 : {TARGET_CHAIN}')
print(f'  Seed residues         : {seed_res_nums}')
print(f'  Sphere radius         : {SPHERE_RADIUS} A')
print(f'  SASA threshold        : {SASA_THRESHOLD_VAL}')
print(f'  Candidate pool size   : {len(candidate_pool)}')
print(f'  Neighbourhood residues: {len(neighbourhood_df)}')
print(f'  Clean PDB             : {clean_pdb_path}')
print()
print('Next: Notebook 00 (Bayesian controller) explores this candidate pool,')
print('and Notebook 02 designs binders for each trial combination using the')
print('clean PDB and hotspot_residues.json.')



# ---- cell f9172326 ----
# Cell 12: Visualise epitope (2D projection)

# Classify residues by amino-acid type for colouring. (Previously computed in
# the old hotspot cell, which has been removed; done here so the plots stand
# alone.)
HYDROPHOBIC_AA = {'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TRP', 'MET', 'TYR'}
CHARGED_AA     = {'ARG', 'LYS', 'ASP', 'GLU', 'HIS'}
POLAR_AA       = {'SER', 'THR', 'ASN', 'GLN', 'CYS'}
if 'aa_type' not in neighbourhood_df.columns:
    neighbourhood_df['aa_type'] = neighbourhood_df['res_name'].apply(
        lambda x: 'hydrophobic' if x in HYDROPHOBIC_AA
        else 'charged'     if x in CHARGED_AA
        else 'polar'       if x in POLAR_AA
        else 'other'
    )

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

COLORS = {'hydrophobic': 'orange', 'charged': 'royalblue', 'polar': 'green', 'other': 'gray'}

# --- Plot 1: XY projection coloured by amino acid type ---
ax = axes[0]
for aa_type, group in neighbourhood_df.groupby('aa_type'):
    ax.scatter(group['x'], group['y'],
               c=COLORS.get(aa_type, 'gray'),
               label=aa_type, s=50, alpha=0.7)

# Highlight core epitope residues
core_rows = neighbourhood_df[neighbourhood_df['in_core']]
if not core_rows.empty:
    ax.scatter(core_rows['x'], core_rows['y'],
               c='red', s=180, marker='*', label='core epitope', zorder=5)
    for _, row in core_rows.iterrows():
        ax.annotate(f"{row['res_name']}{row['res_num']}",
                    (row['x'], row['y']), fontsize=7,
                    xytext=(5, 5), textcoords='offset points')

# Mark epitope COM
ax.scatter(epitope_com[0], epitope_com[1],
           c='black', s=220, marker='+', linewidths=3, label='epitope COM', zorder=6)

ax.set_xlabel('X (Å)')
ax.set_ylabel('Y (Å)')
ax.set_title(f'{TARGET_NAME} — epitope XY projection')
ax.legend(fontsize=8)
ax.set_aspect('equal')

# --- Plot 2: Distance distribution ---
ax = axes[1]
for aa_type, group in neighbourhood_df.groupby('aa_type'):
    ax.scatter(group['dist_to_epitope'],
               range(len(group)),
               c=COLORS.get(aa_type, 'gray'),
               label=aa_type, s=30, alpha=0.7)

ax.axvline(x=8.0,  color='red',    linestyle='--', alpha=0.5, label='hotspot cutoff (8 Å)')
ax.axvline(x=EPITOPE_RADIUS, color='purple', linestyle=':', alpha=0.5,
           label=f'neighbourhood limit ({EPITOPE_RADIUS} Å)')

ax.set_xlabel('Distance from epitope COM (Å)')
ax.set_ylabel('Residue index (sorted)')
ax.set_title('Distance distribution')
ax.legend(fontsize=8)

plt.suptitle(f'{TARGET_NAME} — Epitope Analysis', fontsize=13, fontweight='bold')
plt.tight_layout()

fig_path = output_dir / f'{TARGET_LABEL}_epitope_analysis.png'
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
pass
print(f'Figure saved: {fig_path}')



print('FINAL CANDIDATE POOL:', sorted(candidate_pool))
