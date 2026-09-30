# Basement Antibody Works (BAW)

Open-source computational nanobody design pipeline targeting *Trypanosoma cruzi* trans-sialidase for Chagas disease research. Built and run by a single scientist on consumer hardware. Everything is free and open.

The pipeline goes from a crystal structure of the target protein to scored and ranked nanobody candidate sequences, using RFdiffusion, ProteinMPNN, and RoseTTAFold2. A Bayesian optimization controller searches the space of epitope residue combinations across overnight GPU runs.

Follow the project at [basementantibodyworks.substack.com](https://basementantibodyworks.substack.com)

---

## What is in this repository

```
BasementAntibodyWorks/
  notebooks/     the four pipeline notebooks
  data/ts/       target analysis outputs (hotspot residues, SASA consensus, cleaned 1MS3 structure)
  structures/    figures and structure files
  baw_environment_gtx1070.yml   conda environment (use this one)
  baw_environment.yml           original export from the RTX 4070 laptop (reference only)
```

Trial outputs from the optimization campaign (`designs/`) are not stored in the repository because the files are large. The results table further down lists the scores.

---

## Pipeline

```
Notebook 00a        Notebook 01         Notebook 02                          Notebook 03
(SASA analysis)  →  (Target analysis) → (RFdiffusion + ProteinMPNN + RF2) →  (Bayesian optimization)
```

- **Notebook 00a** and **Notebook 01** select the target surface. Their outputs are already in `data/ts/`, so you can skip them and start at Notebook 03.
- **Notebook 02** runs one design trial: RFdiffusion backbones, ProteinMPNN sequences, RoseTTAFold2 scoring.
- **Notebook 03** decides which residue combinations to test, calls Notebook 02 for each one, and learns from the scores.

---

## Requirements

**Hardware**
- NVIDIA GPU with at least 8 GB of VRAM. The scoring step uses about 95 percent of 8 GB on a GTX 1070, so smaller cards are unlikely to work.
- At least 30 GB of free disk space for the conda environment, RFantibody, and model weights. Allow more for trial outputs.
- Time: the notebook estimates about 3.5 hours per trial on a GTX 1070 with the default settings (10 backbones, 4 sequences each), so a batch of 10 trials takes more than a day. Faster GPUs are quicker.

**Operating system**
- Tested on Ubuntu under WSL2 on Windows, on an RTX 4070 laptop and a GTX 1070 desktop.
- Not tested on native Linux or macOS. Some notebook features assume WSL, see Troubleshooting.

**Software**
- git and curl
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- A current NVIDIA driver. On WSL, install the driver on Windows, not inside Ubuntu. Running `nvidia-smi` inside Ubuntu should show your GPU.

**Important for WSL users:** keep the repository, conda, and RFantibody inside the Linux filesystem (your home folder). Do not install them under `/mnt/c` or `/mnt/d`. Those drives are case-insensitive and conda package extraction fails there.

---

## Installation

### 1. Clone the repository

The notebooks expect the repository at `~/BasementAntibodyWorks`.

```bash
git clone https://github.com/sardism/sialidaseantibodies ~/BasementAntibodyWorks
cd ~/BasementAntibodyWorks
```

### 2. Create the conda environment

Use `baw_environment_gtx1070.yml`. It was exported from a working environment and is the file to use on a fresh machine.

```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
conda env create -f baw_environment_gtx1070.yml
conda activate BAW
```

Skip the two `conda tos` lines if your conda version does not have that command. Creating the environment downloads several GB and takes 10 to 20 minutes.

Check that PyTorch sees the GPU:

```bash
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
```

The second value should be `True`.

### 3. Install RFantibody

RFantibody is not part of the conda environment. It has its own uv-managed environment and must be installed in your home folder:

```bash
cd ~
git clone https://github.com/RosettaCommons/RFantibody.git
cd ~/RFantibody
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.bashrc
bash include/download_weights.sh
uv sync
```

The weights download is several GB. Notebook 02 expects RFantibody at `~/RFantibody` and uv at `~/.local/bin/uv`.

Check the installation:

```bash
cd ~/RFantibody
uv run rfdiffusion --help | head -5
uv run proteinmpnn --help | head -5
uv run rf2 --help | head -5
```

Each command should print a usage message.

### 4. Start JupyterLab

```bash
conda activate BAW
cd ~/BasementAntibodyWorks
jupyter lab
```

Open the link that JupyterLab prints in your browser.

---

## Running the pipeline

### Option A: run the optimization campaign (Notebook 03)

This is the main entry point. It needs `data/ts/ts_hotspot_residues.json` and the cleaned `1MS3_chainA_clean.pdb` next to it, which are included in the repository.

1. **Set the paths.** Open `03_bayesian_controller.ipynb` and run the imports and the configuration cell. The configuration form needs three paths: Notebook 02 (`notebooks/02_backbone_designv2.ipynb`), the designs folder (`designs`, created if missing), and the hotspot file (`data/ts/ts_hotspot_residues.json`).
   - On WSL you can click Confirm and choose the three items in the Windows file dialogs. Select the `designs` folder itself, not a subfolder.
   - On any other system, do not click Confirm. Edit the `NOTEBOOK02_PATH`, `DESIGNS_DIR`, and `HOTSPOT_JSON` variables in the configuration cell instead. The default for `NOTEBOOK02_PATH` points to an older filename, so change it to `02_backbone_designv2.ipynb`.
2. **Load the search space.** The next cell reads the candidate residues and builds all 31,179 valid combinations (3 to 8 residues, with residues 119 and 312 always included). It adds residue contact scores automatically if they are missing.
3. **Load existing trials.** The following cell scans `designs/` for finished `trial_result.json` files. A fresh clone has none.
4. **Cold start.** With fewer than 10 trials, the Latin Hypercube Sampling cell selects 10 evenly spread combinations.
5. **Run the batch.** The run cell shows a table and a Run batch button. Clicking it calls Notebook 02 once per combination and saves results under `designs/trial_XXXX/`. It is safe to interrupt and rerun: finished combinations are skipped.
6. **Learn and repeat.** After at least 3 scored trials, run the Gaussian Process fit, then the acquisition cell (selects the next batch), then the dashboard and uncertainty map, then the run cell again. Repeat until the convergence check says the best score has stopped improving.

Expected output: one folder per trial in `designs/`, each with a `trial_result.json` and the intermediate design files. Results vary from run to run because the design steps are stochastic, so your scores will not match the table below exactly.

### Option B: run a single design trial (Notebook 02)

Open `02_backbone_designv2.ipynb` and run it from the top. In standalone mode it asks for the hotspot file and an output folder, and designs against all candidate residues. The key settings are in its configuration cell:

| Setting | Default | Meaning |
|---|---|---|
| `N_DESIGNS_RFD` | 10 | RFdiffusion backbones per trial |
| `MPNN_SEQS_PER_STRUCT` | 4 | ProteinMPNN sequences per backbone |
| `RF2_RECYCLES` | 10 | RoseTTAFold2 recycling iterations |
| `FILTER_PAE_MAX` | 13.5 | interface pAE cutoff for passing designs |
| `FILTER_RMSD_MAX` | 2.0 | CDR RMSD cutoff in angstroms |

Lower `N_DESIGNS_RFD` or `MPNN_SEQS_PER_STRUCT` to shorten a trial, but very small runs may produce no designs that pass the filters.

### Option C: reproduce the target selection (Notebooks 00a and 01)

`00a_structural_analysis.ipynb` runs the multi-structure surface accessibility analysis and writes `ts_consensus_sasa.json`. `01_target_analysis.ipynb` reads it and writes `ts_hotspot_residues.json` and the cleaned structure. Both write under `data/ts/`. The outputs are already in the repository, so this step is only needed to reproduce the target analysis or to adapt the pipeline to another protein.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `conda env create` complains about Terms of Service | Run the two `conda tos accept` commands from step 2 |
| Conda fails extracting packages with a case-insensitive filesystem error | Move everything into the Linux home folder. Do not use `/mnt/c` or `/mnt/d` |
| pip cannot find `torch==2.7.1+cu118` | You used `baw_environment.yml`. Use `baw_environment_gtx1070.yml` |
| `powershell.exe` not found when clicking Confirm | You are not on WSL. Skip Confirm and set the three path variables by hand |
| File dialog returns paths starting with `//wsl.localhost/...` and loading fails | The path cleanup assumes the WSL distribution is named `Ubuntu`. Edit the `fix_wsl_path` prefixes in the candidate-pool cell to match your distribution name |
| `RFantibody is not installed` or `FileNotFoundError` for `uv` | RFantibody must be at `~/RFantibody` and uv at `~/.local/bin/uv`. Re-run step 3 |
| `FileNotFoundError` for `1MS3_chainA_clean.pdb` | Make sure `data/ts/1MS3_chainA_clean.pdb` exists next to the hotspot JSON. Notebook 02 searches there if the path stored in the JSON does not exist |
| Widget form does not render (`ChunkLoadError`) | Restart JupyterLab and hard refresh the browser. If it persists, run `pip install --upgrade ipywidgets jupyterlab_widgets` and restart |
| CUDA out of memory | Close other GPU programs. RoseTTAFold2 needs close to 8 GB. Cards with less memory are unlikely to work |
| Trials finish with no score | No design passed the filters. Check the RF2 score file in the trial folder and adjust `FILTER_PAE_MAX` |
| Ubuntu will not start and the disk is full | WSL stores its disk image on the drive where it was installed. Free space there or move WSL to a larger drive |

---

## Current results

Interaction pAE from RoseTTAFold2 (lower is better). Trial 1 is the pilot run using all 20 candidate residues. Trials 2 to 11 are the Latin Hypercube Sampling batch.

| Trial | Residues | Best pAE |
|---|---|---|
| 1 (pilot) | all 20 candidate residues | 7.36 |
| 2 | 55, 119, 248, 249, 281, 308, 309, 312 | 7.57 |
| 3 | 55, 118, 119, 248, 249, 308, 312, 313 | **4.15** |
| 4 | 117, 119, 122, 201, 248, 249, 309, 312 | 8.09 |
| 5 | 55, 117, 118, 119, 201, 248, 281, 312 | 12.14 |
| 6 | 13, 117, 118, 119, 122, 249, 311, 312 | 8.93 |
| 7 | 118, 119, 201, 202, 249, 308, 311, 312 | 10.69 |
| 8 | 119, 248, 281, 309, 311, 312, 362 | 11.17 |
| 9 | 118, 119, 202, 248, 281, 308, 309, 312 | 4.65 |
| 10 | 117, 119, 201, 202, 281, 311, 312, 313 | 10.98 |
| 11 | 13, 55, 119, 122, 202, 248, 249, 312 | 5.18 |

- **Target:** *Trypanosoma cruzi* trans-sialidase (PDB: 1MS3)
- **Epitope core identified:** residues 119, 248, 312
- **Campaign status:** LHS batch complete, Gaussian Process guided search in progress
- **Hardware:** Alienware RTX 4070 (primary) and desktop GTX 1070 (secondary)

These are computational predictions. Until validated experimentally they are not molecules.

---

## Links

- Substack: [basementantibodyworks.substack.com](https://basementantibodyworks.substack.com)
- GitHub: [github.com/sardism/sialidaseantibodies](https://github.com/sardism/sialidaseantibodies)

---

## License

MIT. The software is provided as is, without warranty of any kind.

---

## Citation

If you use this pipeline or build on this work please cite the Substack posts and link to this repository. The project is open science. The point is for others to use it.
