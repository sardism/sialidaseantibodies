import json

NB = "/home/sardism/BasementAntibodyWorks/notebooks/02_backbone_designv2.ipynb"
OUT = "/home/sardism/BasementAntibodyWorks/data/dry_run/nb02_dry_run.py"

INJECT_TOP = '''
# === INJECTED (dry run) ===
from pathlib import Path
TRIAL_COMBINATION = [119, 311, 312, 245, 35]
CANDIDATE_POOL_JSON = str(Path("~/BasementAntibodyWorks/data/dry_run/nb01/ts_hotspot_residues.json").expanduser())
# === END INJECTED ===
'''

# N_DESIGNS_RFD etc. must be injected AFTER Cell 2 (Configuration) runs --
# Cell 2 defines its own hardcoded defaults for these names, so injecting
# them earlier (e.g. right after Cell 0) gets silently clobbered when Cell 2
# executes. Found this the hard way: a first dry-run pass showed "RFdiffusion
# N: 50" in the printed config instead of the injected 2.
INJECT_AFTER_CONFIG = '''
# === INJECTED (dry run) ===
output_dir = Path("~/BasementAntibodyWorks/data/dry_run/nb02").expanduser()
output_dir.mkdir(parents=True, exist_ok=True)
N_DESIGNS_RFD=2
MPNN_SEQS_PER_STRUCT=1
RF2_RECYCLES=1
N_RESCUED=2
DETERMINISTIC=True
# === END INJECTED ===
'''

# Extra setup not literally in the instructions but required to make the mock
# self-consistent: there is no real RFantibody install or `uv` binary on this
# machine, so Cell 6's plain filesystem checks (rfab_dir.exists(), .venv,
# weights/*.pt) and Cell 7's framework-file check would fail before any
# subprocess call even happens. We fabricate a minimal fake RFantibody
# directory tree and point RFANTIBODY_DIR at it, then monkeypatch
# subprocess.run / subprocess.Popen so every `uv run <tool> ...` invocation
# is intercepted and produces believable fake Quiver/score/PDB files instead
# of touching the network or a real install.
INJECT_MOCK = r'''
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
'''

data = json.load(open(NB))
cells = data['cells']

def cell_src(c):
    return ''.join(c.get('source', []))

def starts_with(c, prefix):
    return cell_src(c).lstrip().startswith(prefix)

out_parts = []
out_parts.append("# AUTO-GENERATED dry-run script for 02_backbone_designv2.ipynb\n")
out_parts.append("# Popup/controller cells replaced with injected values + mocked subprocess.\n\n")

SKIP_PREFIXES = (
    '# Cell 3: File selection',
)

for c in cells:
    if c.get('cell_type') != 'code':
        continue
    src = cell_src(c)

    if starts_with(c, '# Cell 0: Controller input'):
        out_parts.append(src)
        out_parts.append("\n\n")
        out_parts.append(INJECT_TOP)
        out_parts.append("\n")
        continue

    if any(starts_with(c, p) for p in SKIP_PREFIXES):
        continue

    if starts_with(c, '# Cell 2: Configuration'):
        out_parts.append(src)
        out_parts.append("\n\n")
        out_parts.append(INJECT_AFTER_CONFIG)
        out_parts.append("\n")
        out_parts.append(INJECT_MOCK)
        out_parts.append("\n")
        continue

    out_parts.append(f"\n# ---- {src.splitlines()[0]} ----\n")
    out_parts.append(src)
    out_parts.append("\n\n")

with open(OUT, "w") as f:
    f.write("".join(out_parts))

print("wrote", OUT)
