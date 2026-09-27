Here is the complete README as plain text ready to paste:

---

# Basement Antibody Works (BAW)

Open-source computational nanobody design pipeline targeting *Trypanosoma cruzi* trans-sialidase for Chagas disease research. Built and run by a single scientist on consumer hardware. Everything is free and open.

The pipeline goes from a crystal structure of the target protein all the way to scored and ranked nanobody candidate sequences using RFdiffusion, ProteinMPNN, and RoseTTAFold2, orchestrated by a Bayesian optimization controller that searches the epitope combination space intelligently across overnight GPU runs.

Follow the project at basementantibodyworks.substack.com

---

## Pipeline

Notebook 00a → Notebook 01 → Notebook 02 → Notebook 03

(SASA analysis) → (Target analysis) → (RFdiffusion + ProteinMPNN + RF2) → (Bayesian optimization)

---

## Quick Start

git clone https://github.com/sardism/sialidaseantibodies
cd sialidaseantibodies
conda env create -f baw_environment.yml
conda activate BAW
jupyter lab

Then open the notebooks in order.

---

## Requirements

- Python 3.10
- NVIDIA GPU with CUDA support (tested on RTX 4070 and GTX 1070)
- RFantibody installed separately: https://github.com/RosettaCommons/RFantibody
- Conda environment: baw_environment.yml (RTX 4070 / CUDA 11.8) or baw_environment_gtx1070.yml (GTX 1070 / CUDA 12.6)

---

## Notebooks

00a_structural_analysis.ipynb
Multi-structure Shrake-Rupley SASA analysis across 12 T. cruzi trans-sialidase crystal structures. Mutation filter. Consensus exposure classification. Paratope contact scoring from Reis et al. 2022.

01_target_analysis.ipynb
Geometric sphere search around active site seed residues. Hotspot residue selection with force-include mechanism for biologically essential residues.

02_backbone_designv2.ipynb
RFdiffusion backbone design, ProteinMPNN sequence design, RoseTTAFold2 structural validation and scoring. Logistic regression rescue. Accepts trial combinations from the Bayesian controller.

03_bayesian_controller.ipynb
Gaussian Process surrogate model with Expected Improvement acquisition function. Latin Hypercube Sampling cold start. Live uncertainty map. Resume logic. Papermill integration for unattended overnight runs.

---

## Current Results

Target: Trypanosoma cruzi trans-sialidase (PDB: 1MS3)
Best interaction pAE: 4.15 (trial 3, combination: 55, 118, 119, 248, 249, 308, 312, 313)
Epitope core identified: residues 119, 248, 312
Campaign status: LHS batch complete (11 trials), GP-guided search in progress
Hardware: Alienware RTX 4070 (primary) + desktop GTX 1070 (secondary)

---

## Links

Substack: basementantibodyworks.substack.com
GitHub: github.com/sardism/sialidaseantibodies

---

## License

MIT

---

## Citation

If you use this pipeline or build on this work please cite the Substack posts and link to this repository. The project is open science — the point is for others to use it.
