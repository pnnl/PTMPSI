# PTM-Psi

<p align="center">
<img alt="pnnl logo" src="./.docs/logos/pnnl_logo.png" width="200pt" height="100pt"/> &emsp;
<img alt="pnnl logo" src="./.docs/logos/doe_logo.png" width="200pt" height="60pt"/>
</p>

<br /><br />

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![python](https://img.shields.io/badge/Python-3.9-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![stability-alpha](https://img.shields.io/badge/stability-alpha-f4d03f.svg)](#PTM-Psi)


## Overview

*PTM-Psi* is a Python Package to Facilitate the Computational Investigation of Post-Translational Modification on Protein Structures and Their Impacts on Dynamics and Functions. 

## Installation

### Virtual Environment

The following instructions can be used to install *PTM-Psi* a virtual environment in editable mode. In this way, updates pulled from the repository will become available without the need to reinstall the package.

```bash
python -m venv ptmpsi
source ptmpsi/bin/activate
pip install --upgrade pip

git clone https://github.com/pnnl/PTMPSI.git
cd PTMPSI
pip install -e .
```

### ptmpsi-PCA Analysis Module
This module is designed for performing Principal Component Analysis (PCA) on molecular dynamics (MD) simulations within the PTMPSI framework. It provides a comprehensive workflow covering feature extraction, PCA computation, conformational clustering, scoring, and figure generation. The analysis produces both per-system and combined results, enabling in-depth insights into MD trajectories.

## 🚀 Getting Started
### Navigate to the Module Directory
Begin by changing into the `ptmpsi-PCA` directory:
```
cd PTMPSI/ptmpsi-PCA
```

### Set Up the Environment
Run the setup script to prepare the necessary environment and directories. This step establishes a Conda environment named `gppenv` and only needs to be executed once.
```bash
bash setup.sh  # The "setup.sh" is ran just once to establish the conda environment named "gppenv"
conda activate gppenv 
```

### Generate MD Analysis Pipelines
Execute the driver script to create post-processing directories for your MD systems from the raw simulation data. You will be prompted to provide the following paths and values:
```bash
python generate_pipelines.py 
```
-   `<path_to_RAWsimulations_DIR>`: Path to your raw MD simulation directories.
-   `<numeric_value_of_REFerence_DIR>`: Numeric identifier for the reference directory.
-   `<path_to_POSTprocess_dir>`: Desired path for the generated post-processing directories.

### Run MD Job Analysis
Navigate to the generated workflow directory for your specific system and execute the `md_bundle_flux.sh` script. This script orchestrates the entire MD analysis workflow:
**Navigate to the generated workflow directory** and run the MD jobs:
```bash
cd <generated_workflow_dir>
bash md_bundle_flux.sh
```

#### About `md_bundle_flux.sh`:
This script performs a series of sequential analyses:

-   **Solute Extraction**: Extracts the solute from MD `<.xtc>` trajectory files.
-   **Chain ID Assignment**: Adds chain IDs to the respective reference `<.pdb>` file for each system.
-   **Feature Extraction**: Performs initial feature extraction from the MD data.
-   **PCA & Clustering**: Computes Principal Component Analysis and conducts conformational clustering.
-   **Kinetic Network Generation**: Draws kinetic transition networks based on Mean First Passage Time (MFPT) and cluster residence time.
-   **System Scoring**: Assigns a score to each system based on its PC value. A combined score is also calculated, incorporating MFPT and cluster residence time.
-   **Representative PDBs**: Generates representative PDB files for identified clusters.
-   **Result Generation**: Produces both per-system and combined results, referenced against a combined standard.

### Restarting Analysis (If Required)
If you need to re-run the analysis for any reason, first remove the `.done` flags from the system directories, then re-execute the `md_bundle_flux.sh` script:
**If restart is require**
```bash
find . -maxdepth 2 -type f -name ".done" \( -regex '\./[0-9]+/\.done' -o -regex '\./REF/\.done' \) -delete
```
### Optional: Clean up and Re-run MD Analysis
To clean up previous analysis figures and re-run the MD jobs, follow these steps:

### Check which MD analysis figure directories exist (Optional)

This command lists all `figures` or `figures_*` directories within the system subdirectories (up to two levels deep), helping you verify what will be deleted.
```bash
find . -maxdepth 2 -type d \( -name "figures" -o -name "figures_*" \) \( -regex '\./[0-9]+/figures.*' -o -regex '\./REF/figures.*' \) -print 
```
### Remove MD analysis figures (Optional)
This command deletes all identified `figures` or `figures_*` directories and their contents. Use with caution.
```bash
find . -maxdepth 2 -type d \( -name "figures" -o -name "figures_*" \) \( -regex '\./[0-9]+/figures.*' -o -regex '\./REF/figures.*' \) -exec rm -rf {} +
```
## Re-run MD jobs
After cleaning up, re-execute the main analysis script.
```bash
bash md_bundle_flux.sh
```

---
## 📂 Output
-   **Per-system results**: Located within their respective system directories.
-   **Combined PCA results**:
    -   `PC_score_ranking.csv`
    -   `combined_score_ranking.csv`
-   **Generated Figures**: Visualizations of key MD properties:
    -   RMSD (Root Mean Square Deviation)
    -   Rg (Radius of Gyration)
    -   RMSF (Root Mean Square Fluctuation)
    -   Fraction of Secondary Structure Formation
    -   Contact maps
    -   PCA scatter plots
---

## Citation

## License

*PTM-Psi* is made freely available under the terms of a modified 3-clause BSD license. See [LICENSE](./LICENSE) for details.

## Acknowledgments

The development of *PTM-Psi* is part of the Predictive Phenomics Initiative at the Pacific Northwest National Laboratory (PNNL). Also, a portion of the research was performed using the Molecular Sciences Computing Facility (Tahoma) at the Environmental Molecular Sciences Laboratory (EMSL) and using resources available through Research Computing at PNNL. It was conducted under the Laboratory Directed Research and Development Program at PNNL. PNNL is a multi-program national laboratory operated by Battelle Memorial Institute for the U.S. Department of Energy under contract DE-AC05-76RL01830.
