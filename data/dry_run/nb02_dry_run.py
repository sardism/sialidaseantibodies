# AUTO-GENERATED dry-run script for 02_backbone_designv2.ipynb
# Popup/controller cells replaced with injected values + mocked subprocess.

# Cell 0: Controller input
#
# Notebook 00 (the Bayesian controller) drives this notebook by setting a
# hotspot combination to try and the path to the candidate-pool JSON from
# Notebook 01. For standalone testing, leave both as None and select the JSON
# with the file picker in Cell 3.

import json, sys
from pathlib import Path

# This combination is passed by Notebook 00 (the Bayesian controller)
# For standalone testing, set this manually
TRIAL_COMBINATION = None  # e.g. [119, 245, 311, 312] — set by controller

# Load candidate pool from Notebook 01 output
CANDIDATE_POOL_JSON = None  # path to ts_hotspot_residues.json — set by controller



# === INJECTED (dry run) ===
from pathlib import Path
TRIAL_COMBINATION = [119, 311, 312, 245, 35]
CANDIDATE_POOL_JSON = str(Path("~/BasementAntibodyWorks/data/dry_run/nb01/ts_hotspot_residues.json").expanduser())
# === END INJECTED ===


# ---- # Cell 1: Imports ----
# Cell 1: Imports

import os
import sys
import json
import subprocess
import shutil
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import NeighborSearch

import ipywidgets as widgets
from IPython.display import display, HTML, clear_output

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')

# --------------------------------------------------------------------------
# Popup helpers -- native Windows dialogs via PowerShell (WSL interop)
# --------------------------------------------------------------------------
def ps_browse_file(title='Select file', filter_str='All files|*.*'):
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

# --------------------------------------------------------------------------
# Subprocess helper: run command, stream stdout, raise on failure
# --------------------------------------------------------------------------
def run_cmd(cmd, cwd=None, env=None, label=''):
    """Run a shell command, print output in real time, raise on error."""
    print(f'\n>>> {label}' if label else f'\n>>> Running: {" ".join(str(c) for c in cmd)}')
    print('-' * 60)
    proc = subprocess.Popen(
        [str(c) for c in cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=str(cwd) if cwd else None,
        env=env,
    )
    for line in proc.stdout:
        print(line, end='')
    proc.wait()
    print('-' * 60)
    if proc.returncode != 0:
        raise RuntimeError(f'Command failed with exit code {proc.returncode}')
    return proc.returncode

print('Imports loaded.')


# Cell 2: Configuration
# =============================================================
# Edit these values before running the rest of the notebook.
# =============================================================

# ---  RFantibody repository location  ---
# Where to clone (or find) the RFantibody repo.
RFANTIBODY_DIR = Path.home() / 'RFantibody'

# ---  Antibody format  ---
# 'nanobody' : single heavy chain VHH (default, simpler, more characterized)
# 'scfv'     : single-chain variable fragment (heavy + light chains)
ANTIBODY_FORMAT = 'nanobody'

# ---  CDR loop lengths  ---
# Format: "loop:length" or "loop:min-max"
# Nanobody (H1, H2, H3 only):
CDR_LOOPS_NANOBODY = "H1:7,H2:6,H3:5-13"
# Full scFv (heavy + light):
CDR_LOOPS_SCFV     = "H1:7,H2:6,H3:5-13,L1:8-13,L2:7,L3:9-11"

# ---  Design scale  ---
N_DESIGNS_RFD     = 50    # RFdiffusion backbones (use 50-100 for pilot, 1000+ for production)
MPNN_SEQS_PER_STRUCT = 4  # ProteinMPNN sequences per backbone
RF2_RECYCLES      = 10    # RF2 recycling iterations (10 recommended)

# ---  Target cropping  ---
# Residues within CROP_RADIUS of the epitope COM are kept.
# Strongly recommended: reduces compute and avoids irrelevant structure.
CROP_TARGET       = True
CROP_RADIUS       = 30.0  # Angstroms around the epitope COM

# ---  Filtering thresholds  ---
FILTER_PAE_MAX    = 10.0  # RF2 pAE < this (lower = better)
FILTER_RMSD_MAX   =  2.0  # Design vs RF2-predicted RMSD < this (Angstroms)

# ---  Top designs to export for Notebook 03  ---
N_TOP_DESIGNS     = 10

# ---  Logistic regression rescue (Kothiwal et al. 2026)  ---
# Number of ADDITIONAL designs to rescue beyond the threshold-passing set.
# These are designs that failed one metric but score high on the overall
# feature pattern learned by the LR. Set to 0 to disable.
N_RESCUED         = 5

# ---  Deterministic mode (reproducible results)  ---
DETERMINISTIC     = False  # True for debugging; False for diversity

# Derived: pick the right CDR loop string
CDR_LOOPS = CDR_LOOPS_NANOBODY if ANTIBODY_FORMAT == 'nanobody' else CDR_LOOPS_SCFV

print('Configuration loaded.')
print(f'  Format           : {ANTIBODY_FORMAT}')
print(f'  CDR loops        : {CDR_LOOPS}')
print(f'  RFdiffusion N    : {N_DESIGNS_RFD}')
print(f'  ProteinMPNN seqs : {MPNN_SEQS_PER_STRUCT} per backbone')
print(f'  RF2 recycles     : {RF2_RECYCLES}')
print(f'  Target crop      : {CROP_TARGET} ({CROP_RADIUS} A)')
print(f'  Filter pAE < {FILTER_PAE_MAX}, RMSD < {FILTER_RMSD_MAX}')



# === INJECTED (dry run) ===
output_dir = Path("~/BasementAntibodyWorks/data/dry_run/nb02").expanduser()
output_dir.mkdir(parents=True, exist_ok=True)
N_DESIGNS_RFD=2
MPNN_SEQS_PER_STRUCT=1
RF2_RECYCLES=1
N_RESCUED=2
DETERMINISTIC=True
# === END INJECTED ===


# === INJECTED (dry run) - fake RFantibody install + subprocess mocking ===
import subprocess as _subprocess
import shutil as _shutil

FAKE_RFAB_DIR = output_dir / "fake_rfantibody"
if FAKE_RFAB_DIR.exists():
    _shutil.rmtree(FAKE_RFAB_DIR)
(FAKE_RFAB_DIR / ".venv").mkdir(parents=True, exist_ok=True)
(FAKE_RFAB_DIR / "weights").mkdir(parents=True, exist_ok=True)
(FAKE_RFAB_DIR / "weights" / "dummy.pt").write_text("fake weights\n")
_fw_dir = FAKE_RFAB_DIR / "scripts" / "examples" / "example_inputs"
_fw_dir.mkdir(parents=True, exist_ok=True)
(_fw_dir / "h-NbBCII10.pdb").write_text("REMARK fake nanobody framework\n")
(_fw_dir / "hu-4D5-8_Fv.pdb").write_text("REMARK fake scfv framework\n")

RFANTIBODY_DIR = FAKE_RFAB_DIR

def _extract_flag(args, flag):
    args = list(args)
    if flag in args:
        i = args.index(flag)
        if i + 1 < len(args):
            return args[i + 1]
    return None

def _fake_dispatch(cmd, cwd=None):
    cmd = [str(c) for c in cmd]
    if len(cmd) < 3 or cmd[0] != 'uv' or cmd[1] != 'run':
        return 0, [f'(mock) unrecognized command, treated as no-op: {" ".join(cmd)}']

    tool = cmd[2]
    args = cmd[3:]
    cwd = Path(cwd) if cwd else Path.cwd()

    if tool == 'rfdiffusion' and '--help' in args:
        return 0, ['usage: rfdiffusion [options]  (mocked --help)']

    if tool == 'rfdiffusion':
        out_qv = _extract_flag(args, '--output-quiver')
        n = int(_extract_flag(args, '--num-designs') or N_DESIGNS_RFD)
        tags = [f'design_{i:04d}' for i in range(n)]
        Path(out_qv).write_text('\n'.join(tags) + '\n')
        return 0, [f'(mock) rfdiffusion wrote {n} backbone(s) to {out_qv}']

    if tool == 'proteinmpnn':
        in_qv  = _extract_flag(args, '--input-quiver')
        out_qv = _extract_flag(args, '--output-quiver')
        seqs_per = int(_extract_flag(args, '--seqs-per-struct') or MPNN_SEQS_PER_STRUCT)
        in_tags = Path(in_qv).read_text().strip().splitlines()
        out_tags = [f'{t}_seq{j}' for t in in_tags for j in range(seqs_per)]
        Path(out_qv).write_text('\n'.join(out_tags) + '\n')
        return 0, [f'(mock) proteinmpnn wrote {len(out_tags)} sequence(s) to {out_qv}']

    if tool == 'rf2':
        in_qv  = _extract_flag(args, '--input-quiver')
        out_qv = _extract_flag(args, '--output-quiver')
        in_tags = Path(in_qv).read_text().strip().splitlines()
        Path(out_qv).write_text('\n'.join(in_tags) + '\n')
        return 0, [f'(mock) rf2 predicted {len(in_tags)} structure(s) to {out_qv}']

    if tool == 'qvls':
        qv_path = Path(args[0])
        tags = qv_path.read_text().strip().splitlines() if qv_path.exists() else []
        return 0, tags

    if tool == 'qvscorefile':
        qv_path = Path(args[0])
        tags = qv_path.read_text().strip().splitlines() if qv_path.exists() else []
        sc_path = cwd / f'{qv_path.stem}.sc'
        # Deterministic fake scores: first tag always passes the default
        # filters (pAE < 10, RMSD < 2), the rest fail -- gives the notebook
        # both a "passed filter" design and a "rescue candidate" to exercise.
        lines = ['tag\tpae\trmsd']
        for i, tag in enumerate(tags):
            if i == 0:
                pae, rmsd = 6.0, 1.0
            else:
                pae, rmsd = 15.0 + i, 3.0 + 0.5 * i
            lines.append(f'{tag}\t{pae}\t{rmsd}')
        sc_path.write_text('\n'.join(lines) + '\n')
        return 0, [f'(mock) qvscorefile wrote {sc_path}']

    if tool == 'qvextract':
        qv_path = Path(args[0])
        out_dir_idx = args.index('-o') + 1
        out_dir = Path(args[out_dir_idx])
        out_dir.mkdir(parents=True, exist_ok=True)
        tags = qv_path.read_text().strip().splitlines() if qv_path.exists() else []
        for tag in tags:
            (out_dir / f'{tag}.pdb').write_text(f'REMARK fake design {tag}\n')
        return 0, [f'(mock) qvextract wrote {len(tags)} PDB(s) to {out_dir}']

    return 0, [f'(mock) unrecognized rfab tool "{tool}", treated as no-op']


class _FakeStdout:
    def __init__(self, lines):
        self._lines = [l + '\n' for l in lines]
    def __iter__(self):
        return iter(self._lines)


class _FakePopen:
    def __init__(self, cmd, stdout=None, stderr=None, text=True, cwd=None, env=None):
        rc, lines = _fake_dispatch(cmd, cwd)
        self.returncode = rc
        self.stdout = _FakeStdout(lines)
    def wait(self):
        return self.returncode


class _FakeCompleted:
    def __init__(self, returncode, lines):
        self.returncode = returncode
        self.stdout = '\n'.join(lines) + ('\n' if lines else '')
        self.stderr = ''


def _fake_run(cmd, cwd=None, capture_output=False, text=True, env=None, **kw):
    rc, lines = _fake_dispatch(cmd, cwd)
    return _FakeCompleted(rc, lines)


_subprocess.Popen = _FakePopen
_subprocess.run = _fake_run
# === END INJECTED ===


# ---- # Cell 4: Load Notebook 01 outputs ----
# Cell 4: Load Notebook 01 outputs
#
# The candidate-pool JSON path comes from the controller (CANDIDATE_POOL_JSON,
# Cell 0) if set, otherwise from the file picker in Cell 3.
#
# Hotspot residues to build against:
#   - If TRIAL_COMBINATION is set (by the controller), use it directly.
#   - Otherwise fall back to the candidate pool in the JSON (standalone testing).

json_path = CANDIDATE_POOL_JSON or hotspot_json_path

with open(json_path) as f:
    hotspot_data = json.load(f)

target_name       = hotspot_data['target_name']
source_pdb_name   = hotspot_data['source_pdb']
target_chain      = hotspot_data['chain']
epitope_com       = np.array(hotspot_data['epitope_com'])
clean_pdb_path    = Path(hotspot_data['clean_pdb'])

# Candidate pool from Notebook 01 (new schema); fall back to the old
# 'hotspot_residues' key for backward compatibility.
candidate_pool = hotspot_data.get('candidate_pool',
                                  hotspot_data.get('hotspot_residues', []))

# Decide the hotspot residues for this run
if TRIAL_COMBINATION:
    hotspot_res_nums = list(TRIAL_COMBINATION)
    residue_source   = 'TRIAL_COMBINATION (controller)'
else:
    hotspot_res_nums = list(candidate_pool)
    residue_source   = 'candidate pool from JSON (standalone)'

print(f'NOTEBOOK 01 OUTPUTS LOADED')
print('=' * 60)
print(f'JSON source       : {json_path}')
print(f'Target name       : {target_name}')
print(f'Source PDB        : {source_pdb_name}')
print(f'Chain             : {target_chain}')
print(f'Candidate pool    : {candidate_pool}')
print(f'Residue source    : {residue_source}')
print(f'Hotspot residues  : {hotspot_res_nums}')
print(f'Epitope COM       : ({epitope_com[0]:.2f}, {epitope_com[1]:.2f}, {epitope_com[2]:.2f})')
print(f'Clean PDB         : {clean_pdb_path}')
print()

if not hotspot_res_nums:
    raise ValueError(
        'No hotspot residues to design against. Set TRIAL_COMBINATION in Cell 0 '
        'or ensure the JSON contains a non-empty candidate_pool.'
    )

if not clean_pdb_path.exists():
    raise FileNotFoundError(
        f'Clean PDB not found at {clean_pdb_path}\n'
        f'Run Notebook 01 first, or update the path in the JSON.'
    )

# Build hotspot string for RFantibody: e.g. "A35,A119,A245,A342"
hotspot_str = ','.join(f'{target_chain}{r}' for r in hotspot_res_nums)
print(f'RFantibody hotspot string: {hotspot_str}')



# ---- # Cell 5: Crop target PDB to the epitope region ----
# Cell 5: Crop target PDB to the epitope region
#
# RFdiffusion and RF2 scale as O(N^2) with residue count.
# Cropping to the region around the epitope is strongly recommended.
# Keeps all residues within CROP_RADIUS of the epitope COM.

parser = PDBParser(QUIET=True)
structure = parser.get_structure('target', str(clean_pdb_path))
model     = structure[0]
chain     = model[target_chain]

if CROP_TARGET:
    # Collect protein residues that have a CA atom (skip any with missing/disordered CA)
    ca_residues = [r for r in chain if r.get_id()[0] == ' ' and 'CA' in r]
    ca_coords   = np.array([r['CA'].get_vector().get_array() for r in ca_residues])

    # Find residues within CROP_RADIUS of epitope COM
    dists_to_com = np.linalg.norm(ca_coords - epitope_com, axis=1)
    keep_mask    = dists_to_com <= CROP_RADIUS

    protein_res  = [r for r in chain if r.get_id()[0] == ' ']
    keep_res_ids = set(
        ca_residues[i].get_id()[1]
        for i, keep in enumerate(keep_mask)
        if keep
    )

    print(f'TARGET CROPPING')
    print('=' * 60)
    print(f'Original residues : {len(protein_res)}')
    print(f'Residues within {CROP_RADIUS} A of epitope COM: {len(keep_res_ids)}')
    print(f'Residue range kept: {min(keep_res_ids)} — {max(keep_res_ids)}')
    
    # Verify all hotspot residues are included
    missing_hotspots = [r for r in hotspot_res_nums if r not in keep_res_ids]
    if missing_hotspots:
        print(f'WARNING: hotspot residues not in crop: {missing_hotspots}')
        print(f'Increase CROP_RADIUS in Cell 2 to include them.')
    else:
        print(f'All {len(hotspot_res_nums)} hotspot residues are within the crop.')

    class CropSelect(Select):
        def accept_chain(self, c):
            return c.get_id() == target_chain
        def accept_residue(self, r):
            return r.get_id()[0] == ' ' and r.get_id()[1] in keep_res_ids
        def accept_atom(self, a):
            altloc = a.get_altloc()
            return altloc == ' ' or altloc == 'A'

    target_for_design = output_dir / f'{clean_pdb_path.stem}_cropped.pdb'
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(target_for_design), CropSelect())
    print(f'\nCropped target saved: {target_for_design}')

else:
    target_for_design = clean_pdb_path
    protein_res = [r for r in chain if r.get_id()[0] == ' ']
    print(f'Cropping disabled. Using full target: {len(protein_res)} residues.')
    print(f'NOTE: large targets slow down RFdiffusion significantly.')



# ---- # Cell 6: Check RFantibody installation (preflight only) ----
# Cell 6: Check RFantibody installation (preflight only)
#
# RFantibody uses its own uv-managed Python environment (.venv) and is
# installed separately from the BAW conda environment -- this cell only
# verifies the install is present; it does not clone or download anything.

rfab_dir = Path(RFANTIBODY_DIR)
print('RFANTIBODY PREFLIGHT CHECK')
print('=' * 60)
checks = [
    (rfab_dir.exists(), f'Repo found at {rfab_dir}'),
    ((rfab_dir / '.venv').exists(), 'Python environment (.venv) exists'),
    ((rfab_dir / 'weights').exists() and any((rfab_dir/'weights').glob('*.pt')), 'Model weights present'),
]
all_ok = True
for ok, msg in checks:
    status = 'OK' if ok else 'MISSING'
    print(f'  [{status}] {msg}')
    if not ok: all_ok = False
if not all_ok:
    raise RuntimeError(
        'RFantibody is not installed. Run the one-time install:\n'
        '  cd ~ && git clone https://github.com/RosettaCommons/RFantibody.git\n'
        '  cd ~/RFantibody && curl -LsSf https://astral.sh/uv/install.sh | sh\n'
        '  source ~/.bashrc && bash include/download_weights.sh && uv sync'
    )
result = subprocess.run(['uv','run','rfdiffusion','--help'], cwd=str(rfab_dir), capture_output=True, text=True)
if result.returncode != 0:
    raise RuntimeError('rfdiffusion command failed. Check RFantibody installation.')
print('  [OK] rfdiffusion command verified')
print()
print('RFantibody ready.')

def rfab_cmd(tool, *args):
    return ['uv', 'run', tool] + list(args)



# ---- # Cell 7: Select antibody framework ----
# Cell 7: Select antibody framework
#
# RFantibody requires an HLT-formatted framework PDB as the antibody scaffold.
# The RFantibody repo ships with two ready-to-use frameworks:
#   Nanobody : RFantibody/scripts/examples/example_inputs/h-NbBCII10.pdb
#   ScFv     : RFantibody/scripts/examples/example_inputs/hu-4D5-8_Fv.pdb
#
# You can also download additional frameworks from SAbDab (opig.stats.ox.ac.uk)
# and convert them with: python scripts/util/chothia_to_HLT.py -inpdb in.pdb -outpdb out.pdb

BUNDLED_FRAMEWORKS = {
    'nanobody': rfab_dir / 'scripts' / 'examples' / 'example_inputs' / 'h-NbBCII10.pdb',
    'scfv':     rfab_dir / 'scripts' / 'examples' / 'example_inputs' / 'hu-4D5-8_Fv.pdb',
}

default_framework = BUNDLED_FRAMEWORKS.get(ANTIBODY_FORMAT)

print(f'FRAMEWORK SELECTION')
print('=' * 60)
print(f'Format: {ANTIBODY_FORMAT}')
print(f'Default framework: {default_framework}')

use_default = True  # Set to False and run ps_browse_file() below to use a custom framework

if use_default:
    framework_pdb = default_framework
    if not framework_pdb.exists():
        raise FileNotFoundError(
            f'Bundled framework not found at {framework_pdb}\n'
            f'Verify RFantibody is installed correctly (Cell 6).'
        )
    print(f'Using bundled framework: {framework_pdb}')
else:
    # Uncomment to browse for a custom framework:
    # framework_pdb = Path(ps_browse_file(
    #     title='Select HLT-formatted antibody framework PDB',
    #     filter_str='PDB files|*.pdb|All files|*.*'
    # ))
    pass

print(f'\nFramework: {framework_pdb.name}')
print(f'CDR loops : {CDR_LOOPS}')



# ---- # Cell 8: Run RFdiffusion — backbone design ----
# Cell 8: Run RFdiffusion — backbone design
#
# Generates N_DESIGNS_RFD antibody backbone structures docked to the target,
# with CDR loops positioned to contact the hotspot residues.
#
# Output: 1_rfdiffusion.qv (Quiver file containing all backbone PDBs)
#
# Runtime: ~2-5 min for 10 designs on RTX 4070.
# Checkpoint: skips if 1_rfdiffusion.qv already exists.

rfd_output_qv = output_dir / '1_rfdiffusion.qv'

if rfd_output_qv.exists():
    print(f'Checkpoint: {rfd_output_qv.name} already exists. Skipping RFdiffusion.')
    print(f'Delete {rfd_output_qv} to re-run.')
else:
    cmd = rfab_cmd(
        'rfdiffusion',
        '--target',        str(target_for_design),
        '--framework',     str(framework_pdb),
        '--output-quiver', str(rfd_output_qv),
        '--num-designs',   str(N_DESIGNS_RFD),
        '--design-loops',  CDR_LOOPS,
        '--hotspots',      hotspot_str,
    )
    if DETERMINISTIC:
        cmd.append('--deterministic')

    print(f'RFDIFFUSION — BACKBONE DESIGN')
    print('=' * 60)
    print(f'Target     : {target_for_design.name}')
    print(f'Framework  : {framework_pdb.name}')
    print(f'Hotspots   : {hotspot_str}')
    print(f'CDR loops  : {CDR_LOOPS}')
    print(f'N designs  : {N_DESIGNS_RFD}')
    print(f'Output     : {rfd_output_qv}')
    print()

    run_cmd(cmd, cwd=rfab_dir, label='rfdiffusion backbone design')

# Verify output and report count
result = subprocess.run(
    rfab_cmd('qvls', str(rfd_output_qv)),
    cwd=str(rfab_dir),
    capture_output=True, text=True
)
n_backbones = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
print(f'\nBackbones in Quiver: {n_backbones}')



# ---- # Cell 9: Run ProteinMPNN — CDR sequence design ----
# Cell 9: Run ProteinMPNN — CDR sequence design
#
# Takes the RFdiffusion backbone Quiver and designs amino acid sequences
# for all CDR loops. Generates MPNN_SEQS_PER_STRUCT sequences per backbone.
# Total sequences = N_DESIGNS_RFD × MPNN_SEQS_PER_STRUCT
#
# Output: 2_proteinmpnn.qv
#
# Checkpoint: skips if output already exists.

mpnn_output_qv = output_dir / '2_proteinmpnn.qv'

if mpnn_output_qv.exists():
    print(f'Checkpoint: {mpnn_output_qv.name} already exists. Skipping ProteinMPNN.')
    print(f'Delete {mpnn_output_qv} to re-run.')
else:
    cmd = rfab_cmd(
        'proteinmpnn',
        '--input-quiver',   str(rfd_output_qv),
        '--output-quiver',  str(mpnn_output_qv),
        '--seqs-per-struct', str(MPNN_SEQS_PER_STRUCT),
        '--temperature',    '0.2',
    )
    if DETERMINISTIC:
        cmd.append('--deterministic')

    print('PROTEINMPNN — CDR SEQUENCE DESIGN')
    print('=' * 60)
    print(f'Input        : {rfd_output_qv.name}')
    print(f'Output       : {mpnn_output_qv.name}')
    print(f'Seqs/struct  : {MPNN_SEQS_PER_STRUCT}')
    print(f'Temperature  : 0.2')
    expected_total = n_backbones * MPNN_SEQS_PER_STRUCT
    print(f'Expected total sequences: {expected_total}')
    print()

    run_cmd(cmd, cwd=rfab_dir, label='proteinmpnn sequence design')

result = subprocess.run(
    rfab_cmd('qvls', str(mpnn_output_qv)),
    cwd=str(rfab_dir),
    capture_output=True, text=True
)
n_sequences = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
print(f'\nSequences in Quiver: {n_sequences}')



# ---- # Cell 10: Run RF2 — structure prediction and scoring ----
# Cell 10: Run RF2 — structure prediction and scoring
#
# RF2 (antibody-finetuned RoseTTAFold2) predicts the structure of each
# designed sequence and scores confidence. Key output metrics:
#   pAE   — predicted aligned error at the antibody-target interface (lower = better)
#   RMSD  — backbone RMSD between the designed structure and the RF2 prediction
#
# Recommended filter: pAE < 10 AND RMSD < 2 Å
#
# Runtime: ~1-2 min per sequence on RTX 4070.
# For n_sequences = 200: ~3-6 hours. Plan as an overnight run.
# Checkpoint: skips if output already exists.

rf2_output_qv = output_dir / '3_rf2.qv'

if rf2_output_qv.exists():
    print(f'Checkpoint: {rf2_output_qv.name} already exists. Skipping RF2.')
    print(f'Delete {rf2_output_qv} to re-run.')
else:
    cmd = rfab_cmd(
        'rf2',
        '--input-quiver',   str(mpnn_output_qv),
        '--output-quiver',  str(rf2_output_qv),
        '--num-recycles',   str(RF2_RECYCLES),
    )

    print('RF2 — STRUCTURE PREDICTION AND SCORING')
    print('=' * 60)
    print(f'Input          : {mpnn_output_qv.name}')
    print(f'Output         : {rf2_output_qv.name}')
    print(f'Recycles       : {RF2_RECYCLES}')
    print(f'Total to score : {n_sequences}')
    print()
    print('NOTE: this is the slow step. Budget ~1-2 min per sequence on RTX 4070.')
    mins_estimate = n_sequences * 1.5
    print(f'Estimated runtime: {mins_estimate:.0f} min ({mins_estimate/60:.1f} hours)')
    print('Tip: run this as an overnight job.')
    print()

    run_cmd(cmd, cwd=rfab_dir, label='rf2 structure prediction')

result = subprocess.run(
    rfab_cmd('qvls', str(rf2_output_qv)),
    cwd=str(rfab_dir),
    capture_output=True, text=True
)
n_predicted = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
print(f'\nPredicted structures in Quiver: {n_predicted}')



# ---- # Cell 11: Parse RF2 scores and filter designs ----
# Cell 11: Parse RF2 scores and filter designs
#
# Extracts the score table from the RF2 Quiver and applies filtering thresholds.
# Also extracts the RFdiffusion hotspot recovery scores.

rf2_scorefile  = output_dir / '3_rf2.sc'
rfd_scorefile  = output_dir / '1_rfdiffusion.sc'

# Extract RF2 scores
run_cmd(
    rfab_cmd('qvscorefile', str(rf2_output_qv)),
    cwd=rfab_dir,
    label='extract RF2 scores'
)
# qvscorefile writes to <quiver_stem>.sc by default
default_rf2_sc = rfab_dir / '3_rf2.sc'
if default_rf2_sc.exists() and not rf2_scorefile.exists():
    shutil.copy(default_rf2_sc, rf2_scorefile)

# Extract RFdiffusion hotspot recovery scores
run_cmd(
    rfab_cmd('qvscorefile', str(rfd_output_qv)),
    cwd=rfab_dir,
    label='extract RFdiffusion scores'
)
default_rfd_sc = rfab_dir / '1_rfdiffusion.sc'
if default_rfd_sc.exists() and not rfd_scorefile.exists():
    shutil.copy(default_rfd_sc, rfd_scorefile)

# ---- Load and inspect RF2 scores ----
if rf2_scorefile.exists():
    scores_df = pd.read_csv(rf2_scorefile, sep='\t')
    print(f'RF2 SCORES')
    print('=' * 60)
    print(f'Total designs scored: {len(scores_df)}')
    print(f'Columns: {list(scores_df.columns)}')
    print()
    display(scores_df.describe())
else:
    # Try to find any .sc file produced
    sc_files = list(output_dir.glob('*.sc')) + list(rfab_dir.glob('3_rf2*.sc'))
    if sc_files:
        scores_df = pd.read_csv(sc_files[0], sep='\t')
        rf2_scorefile = sc_files[0]
        print(f'Loaded scores from: {sc_files[0]}')
        display(scores_df.head())
    else:
        print('WARNING: No score file found. Check the RF2 run output above.')
        scores_df = pd.DataFrame()



# ---- # Cell 12: Apply filters and rank designs ----
# Cell 12: Apply filters and rank designs

if scores_df.empty:
    print('No scores to filter. Check Cell 11.')
else:
    # Identify pAE and RMSD columns
    # Column names can vary by RFantibody version; find them robustly
    pae_col  = next((c for c in scores_df.columns if 'pae'  in c.lower()), None)
    rmsd_col = next((c for c in scores_df.columns if 'rmsd' in c.lower()), None)
    tag_col  = next((c for c in scores_df.columns if 'tag'  in c.lower() or 
                     'description' in c.lower() or 'name' in c.lower()), None)

    print(f'Score columns identified:')
    print(f'  pAE column   : {pae_col}')
    print(f'  RMSD column  : {rmsd_col}')
    print(f'  Tag column   : {tag_col}')
    print()

    # Apply filters
    filtered = scores_df.copy()
    if pae_col:
        filtered = filtered[filtered[pae_col] < FILTER_PAE_MAX]
    if rmsd_col:
        filtered = filtered[filtered[rmsd_col] < FILTER_RMSD_MAX]

    print(f'FILTERING RESULTS')
    print('=' * 60)
    print(f'Before filter : {len(scores_df)} designs')
    print(f'After filter  : {len(filtered)} designs')
    print(f'  pAE < {FILTER_PAE_MAX}  AND  RMSD < {FILTER_RMSD_MAX} A')
    print()

    # Rank by pAE (lower is better)
    if pae_col and len(filtered) > 0:
        filtered = filtered.sort_values(pae_col).reset_index(drop=True)
        filtered['rank'] = filtered.index + 1

        print(f'TOP {min(N_TOP_DESIGNS, len(filtered))} DESIGNS BY pAE:')
        display_cols = ['rank']
        if tag_col:  display_cols.append(tag_col)
        if pae_col:  display_cols.append(pae_col)
        if rmsd_col: display_cols.append(rmsd_col)
        display(filtered[display_cols].head(N_TOP_DESIGNS))
    elif len(filtered) == 0:
        print('No designs passed the filter.')
        print('Consider relaxing FILTER_PAE_MAX or FILTER_RMSD_MAX in Cell 2.')
        # `filtered` stays empty on purpose: these designs did NOT pass the
        # quality thresholds and must not be exported as "top designs" or
        # counted as passing in Cell 13 / the summary JSON.
        preview_df = scores_df.copy()
        if pae_col:
            preview_df = preview_df.sort_values(pae_col)
        print(f'Showing top 5 unfiltered designs for inspection (none passed filter):')
        display(preview_df.head(5))

    # Save filtered scores
    filtered_csv = output_dir / 'rf2_scores_filtered.csv'
    filtered.to_csv(filtered_csv, index=False)
    print(f'\nFiltered scores saved: {filtered_csv}')



# ---- # Cell 12b: Logistic regression rescue ----
# Cell 12b: Logistic regression rescue
#
# Based on: Kothiwal et al. (2026) Cell Systems 17, 101645
# Key finding from that paper: there is no clear correlation between final-round
# enrichment score and actual binding potency. The best antibodies came from
# lower-abundance clusters that the hard threshold discarded.
#
# In BAW the equivalent mistake is cutting solely on pAE. A design with pAE = 11
# that just missed the cutoff might have strong hotspot contacts, good BSA, and
# a high ProteinMPNN log-likelihood. This cell trains an L1 logistic regression
# on ALL available scores to identify those near-miss candidates.
#
# Output:
#   rescued_df              — up to N_RESCUED additional designs
#   rf2_scores_lr_scored.csv — every design with its LR score attached

import warnings

rescued_df = pd.DataFrame()

if scores_df.empty:
    print('No scores available. Skipping logistic regression rescue.')
else:
    print('LOGISTIC REGRESSION RESCUE (Kothiwal et al. 2026)')
    print('=' * 60)

    # ── Merge RFdiffusion scores if available ──────────────────────────────
    merged_df = scores_df.copy()
    if rfd_scorefile.exists():
        try:
            rfd_df   = pd.read_csv(rfd_scorefile, sep='\t')
            rfd_tag  = next((c for c in rfd_df.columns
                             if any(k in c.lower() for k in ['tag','name','description'])), None)
            if rfd_tag and tag_col:
                rfd_numeric = rfd_df.select_dtypes(include=[float, int]).columns.tolist()
                rfd_rename  = {c: f'rfd_{c}' for c in rfd_numeric}
                rfd_df      = rfd_df.rename(columns=rfd_rename)
                rfd_df      = rfd_df.rename(columns={rfd_tag: tag_col})
                merged_df   = pd.merge(merged_df, rfd_df[[tag_col] + list(rfd_rename.values())],
                                       on=tag_col, how='left')
                new_cols = list(rfd_rename.values())
                print(f'Merged RFdiffusion scores: {new_cols}')
        except Exception as e:
            print(f'Could not merge RFdiffusion scores: {e}')

    # ── Build feature matrix ───────────────────────────────────────────────
    feature_cols = merged_df.select_dtypes(include=[float, int]).columns.tolist()
    feature_cols = [c for c in feature_cols if merged_df[c].notna().any()]
    print(f'Features: {feature_cols}')
    print(f'Total designs: {len(merged_df)}')

    # ── Label: passed both hard thresholds = positive ─────────────────────
    pos_mask = pd.Series(True, index=merged_df.index)
    if pae_col  and pae_col  in merged_df.columns:
        pos_mask &= (merged_df[pae_col]  < FILTER_PAE_MAX)
    if rmsd_col and rmsd_col in merged_df.columns:
        pos_mask &= (merged_df[rmsd_col] < FILTER_RMSD_MAX)

    n_pos = int(pos_mask.sum())
    n_neg = int((~pos_mask).sum())
    print(f'Positive (passed threshold): {n_pos}')
    print(f'Negative (failed threshold): {n_neg}')

    if n_pos < 2 or n_neg < 2 or not feature_cols:
        print()
        print('Not enough data to train logistic regression.')
        print(f'Need at least 2 positive and 2 negative examples with numeric scores.')
        print('Run with more designs (N_DESIGNS_RFD) to enable this cell.')
    else:
        X = merged_df[feature_cols].fillna(merged_df[feature_cols].median())
        y = pos_mask.astype(int)

        scaler  = StandardScaler()
        X_sc    = scaler.fit_transform(X)

        clf = LogisticRegression(
            penalty='l1', C=1.0,
            class_weight='balanced',
            solver='liblinear',
            random_state=42, max_iter=1000
        )
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            clf.fit(X_sc, y)

        lr_scores          = clf.predict_proba(X_sc)[:, 1]
        merged_df['lr_score'] = lr_scores

        try:
            auc = roc_auc_score(y, lr_scores)
            print(f'\nTraining AUC: {auc:.3f}')
            if auc > 0.7:
                print('  Good signal — the LR has learned a meaningful feature pattern.')
            elif auc > 0.55:
                print('  Moderate signal — treat rescued designs with extra scepticism.')
            else:
                print('  Weak signal — scores may not be informative enough for rescue.')
        except Exception:
            pass

        coef_s   = pd.Series(clf.coef_[0], index=feature_cols)
        nonzero  = (coef_s != 0).sum()
        print(f'Non-zero coefficients: {nonzero} / {len(feature_cols)}')
        if nonzero > 0:
            print('Top positive features:')
            print(coef_s[coef_s > 0].sort_values(ascending=False)
                    .head(5).to_string())

        # ── Select rescued designs ─────────────────────────────────────────
        candidates = merged_df[~pos_mask].copy().sort_values('lr_score', ascending=False)

        if N_RESCUED > 0 and len(candidates) > 0:
            rescued_df = candidates.head(N_RESCUED).copy()
            rescued_df['rescue_method'] = 'lr_rescue'

            print(f'\nRESCUED DESIGNS — top {min(N_RESCUED, len(candidates))} by LR score')
            print('-' * 60)
            show = [c for c in [tag_col, 'lr_score', pae_col, rmsd_col] if c and c in rescued_df.columns]
            display(rescued_df[show].reset_index(drop=True))
        else:
            print('\nNo candidates available for rescue or N_RESCUED = 0.')

        lr_csv = output_dir / 'rf2_scores_lr_scored.csv'
        merged_df.to_csv(lr_csv, index=False)
        print(f'\nFull LR-scored table: {lr_csv}')



# ---- # Cell 13: Export top designs for Notebook 03 ----
# Cell 13: Export top designs for Notebook 03
#
# Exports two sets:
#   top_designs/      — designs that passed the hard pAE + RMSD thresholds
#   rescued_designs/  — additional designs rescued by the LR (Cell 12b)
# Notebook 03 receives both sets and can process them independently.

top_designs_dir     = output_dir / 'top_designs'
rescued_designs_dir = output_dir / 'rescued_designs'
top_designs_dir.mkdir(exist_ok=True)
rescued_designs_dir.mkdir(exist_ok=True)

n_top_exported     = 0
n_rescued_exported = 0

if not scores_df.empty and tag_col and len(filtered) > 0:
    top_n    = min(N_TOP_DESIGNS, len(filtered))
    top_tags = filtered.head(top_n)[tag_col].tolist() if tag_col else []

    print(f'EXPORTING TOP {top_n} THRESHOLD-PASSING DESIGNS')
    print('=' * 60)

    if top_tags:
        extract_cmd = rfab_cmd('qvextract', str(rf2_output_qv), '-o', str(top_designs_dir))
        run_cmd(extract_cmd, cwd=rfab_dir, label='extract top design PDBs')
        top_tag_set = set(top_tags)
        for pdb_file in top_designs_dir.glob('*.pdb'):
            if pdb_file.stem not in top_tag_set:
                pdb_file.unlink()
        n_top_exported = len(list(top_designs_dir.glob('*.pdb')))
        print(f'Extracted {n_top_exported} PDB files to {top_designs_dir}')

# ── Export rescued designs ────────────────────────────────────────────────────
if not rescued_df.empty and tag_col and tag_col in rescued_df.columns:
    rescued_tags = rescued_df[tag_col].tolist()
    print(f'\nEXPORTING {len(rescued_tags)} RESCUED DESIGNS')
    print('=' * 60)

    extract_cmd = rfab_cmd('qvextract', str(rf2_output_qv), '-o', str(rescued_designs_dir))
    run_cmd(extract_cmd, cwd=rfab_dir, label='extract rescued design PDBs')
    rescued_tag_set = set(rescued_tags)
    for pdb_file in rescued_designs_dir.glob('*.pdb'):
        if pdb_file.stem not in rescued_tag_set:
            pdb_file.unlink()
    n_rescued_exported = len(list(rescued_designs_dir.glob('*.pdb')))
    print(f'Extracted {n_rescued_exported} PDB files to {rescued_designs_dir}')
else:
    print('\nNo rescued designs to export.')

# ── Save design summary JSON for Notebook 03 ─────────────────────────────────
summary = {
    'target_name':          target_name,
    'target_chain':         target_chain,
    'clean_pdb':            str(clean_pdb_path),
    'cropped_pdb':          str(target_for_design),
    'hotspot_residues':     hotspot_res_nums,
    'hotspot_str':          hotspot_str,
    'framework':            str(framework_pdb),
    'antibody_format':      ANTIBODY_FORMAT,
    'cdr_loops':            CDR_LOOPS,
    'n_designs_rfd':        N_DESIGNS_RFD,
    'n_sequences':          n_sequences,
    'n_passed_filter':      len(filtered),
    'n_top_exported':       n_top_exported,
    'n_rescued_exported':   n_rescued_exported,
    'filter_pae_max':       FILTER_PAE_MAX,
    'filter_rmsd_max':      FILTER_RMSD_MAX,
    'rf2_scorefile':        str(output_dir / 'rf2_scores_filtered.csv'),
    'lr_scorefile':         str(output_dir / 'rf2_scores_lr_scored.csv'),
    'top_designs_dir':      str(top_designs_dir),
    'rescued_designs_dir':  str(rescued_designs_dir),
    'quiver_rf2':           str(rf2_output_qv),
}

summary_json = output_dir / 'nb02_design_summary.json'
with open(summary_json, 'w') as f:
    json.dump(summary, f, indent=2)
print(f'\nDesign summary JSON: {summary_json}')

print()
print('=' * 60)
print(f'NOTEBOOK 02 COMPLETE — {target_name}')
print('=' * 60)
print()
print(f'  RFdiffusion backbones       : {n_backbones}')
print(f'  ProteinMPNN sequences       : {n_sequences}')
print(f'  RF2 predictions             : {n_predicted}')
print(f'  Passed hard filter          : {len(filtered) if not scores_df.empty else "N/A"}')
print(f'  Top designs exported        : {n_top_exported}')
print(f'  Rescued by LR               : {n_rescued_exported}')
print(f'  Total forward to Notebook 03: {n_top_exported + n_rescued_exported}')
print()
print('Next: Notebook 03 uses nb02_design_summary.json for sequence')
print('optimisation and humanisation with ProteinMPNN.')



# ---- # Cell 14: Visualize score distributions ----
# Cell 14: Visualize score distributions

if not scores_df.empty and (pae_col or rmsd_col):
    fig, axes = plt.subplots(1, 2 if (pae_col and rmsd_col) else 1, figsize=(12, 5))
    if not isinstance(axes, np.ndarray):
        axes = [axes]

    ax_idx = 0

    if pae_col:
        ax = axes[ax_idx]; ax_idx += 1
        ax.hist(scores_df[pae_col].dropna(), bins=30, color='steelblue', alpha=0.7, label='all')
        if len(filtered) > 0:
            ax.hist(filtered[pae_col].dropna(), bins=30, color='orange', alpha=0.7, label='passed filter')
        ax.axvline(FILTER_PAE_MAX, color='red', linestyle='--', label=f'cutoff ({FILTER_PAE_MAX})')
        ax.set_xlabel('RF2 pAE (interface)')
        ax.set_ylabel('Count')
        ax.set_title('pAE distribution')
        ax.legend()

    if rmsd_col and ax_idx < len(axes):
        ax = axes[ax_idx]
        ax.hist(scores_df[rmsd_col].dropna(), bins=30, color='steelblue', alpha=0.7, label='all')
        if len(filtered) > 0:
            ax.hist(filtered[rmsd_col].dropna(), bins=30, color='orange', alpha=0.7, label='passed filter')
        ax.axvline(FILTER_RMSD_MAX, color='red', linestyle='--', label=f'cutoff ({FILTER_RMSD_MAX} Å)')
        ax.set_xlabel('RMSD design vs RF2 prediction (Å)')
        ax.set_ylabel('Count')
        ax.set_title('RMSD distribution')
        ax.legend()

    plt.suptitle(f'{target_name} — RFantibody Design Scores', fontsize=13, fontweight='bold')
    plt.tight_layout()

    fig_path = output_dir / 'nb02_score_distributions.png'
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    pass
    print(f'Figure saved: {fig_path}')

    # pAE vs RMSD scatter (if both available)
    if pae_col and rmsd_col and len(scores_df) > 1:
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(scores_df[rmsd_col], scores_df[pae_col],
                   c='steelblue', alpha=0.5, s=20, label='all designs')
        if len(filtered) > 0:
            ax.scatter(filtered[rmsd_col], filtered[pae_col],
                       c='orange', alpha=0.8, s=40, label='passed filter', zorder=5)
        ax.axvline(FILTER_RMSD_MAX, color='red', linestyle='--', alpha=0.5)
        ax.axhline(FILTER_PAE_MAX,  color='red', linestyle='--', alpha=0.5)
        ax.set_xlabel('RMSD (Å)')
        ax.set_ylabel('pAE')
        ax.set_title(f'{target_name} — pAE vs RMSD')
        ax.legend()
        plt.tight_layout()

        scatter_path = output_dir / 'nb02_pae_vs_rmsd.png'
        plt.savefig(scatter_path, dpi=150, bbox_inches='tight')
        pass
        print(f'Scatter plot saved: {scatter_path}')
else:
    print('No scores available to plot.')



# ---- # Cell 15: Write trial result for the controller ----
# Cell 15: Write trial result for the controller
#
# Notebook 00 (the Bayesian controller) reads trial_result.json to score this
# combination and decide what to try next.

trial_result = {
    'combination': TRIAL_COMBINATION if TRIAL_COMBINATION else hotspot_res_nums,
    'best_pae': float(filtered[pae_col].min()) if (not scores_df.empty and pae_col and len(filtered) > 0) else None,
    'n_passed': len(filtered) if not scores_df.empty else 0,
    'output_dir': str(output_dir),
}
result_path = output_dir / 'trial_result.json'
with open(result_path, 'w') as f:
    json.dump(trial_result, f, indent=2)
print(f'Trial result saved: {result_path}')


