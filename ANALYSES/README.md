## General Information ##
This repository contains the analysis scripts and data released to accompany: Cyanophage CP12 Rewires Host Carbon Regulation through Interface Remodeling and Redox Buffering.

## Project Overview ##
This work was conducted under the NW-BRaVE (Northwest Biopreparedness Research Virtual Environment) project by Pacific Northwest National Laboratory. Our vision is to create a powerful and user-friendly platform to elucidate the fundamental principles of how molecular interactions drive pathogen-host relationships and host shifts. Our overall focus is on elucidating general principles of detecting, assessing, and surveilling molecular interaction, adaptation, and coevolution that are system agnostic and extensible to other viral-host interactions. 

## Installation ##
First clone this repository and then in an environment with ```Python >= 3.8``` install the relevant python packages:
```pip install -r requirements.txt```

Some subfolders document additional requirements of their own. Install those separately when you use that workflow.


## Descriptions and Inventory ##

### CP12 protein family build ###
These scripts compiles a protein family collection for CP12 pulled from meta-genome derived sequences and sequences in public reference databases (NCBI and UniProt). 
* **Data** is located in ```CP12_protein_family_build/data```. **Data Type**: sequences in .tsv, .fasta, .hmm, .txt and .sto

### Static interface HMM mapping ###
These scripts map the protein-partner (alternate complex subunits) static interfaces on to the CP12 family HMM columns. 
* **Data** is located in ```Static_interface_HMM_mapping/data```. Data Type: sequences (.tsv), sequence alignment (.hmmm), and structures (PDB)

### CP12 protein family HMM column metrics ###
These scripts calculate a variety of metrics (sequence diversity, acidity, Shannon entropy, etc.) per column positions across the protein family sequences. 
* **Data** is located in ```CP12_Protein_Family_HMM_Column_Metrics/data```. Data Types include sequence alignments (.hmmm), as well as class and feature predicitons often in .tsv format

### Physico-Chemical Changes ###
Calculates the change in change in volume ($Å^{3}$) hydrophobicity (Kyte-Doolittle), and the charge change (charge at pH 7.4), given the direct sequence alignment between MED4 and PHM2 CP12 sequences. 
* **Data** located in ```Physico-Chemical/data/HMM_seq_align.fasta``` (reduced version of ```cp12_combined_residue_to_hmmcol_mappingOUT.txt``` with two entries of just MED4 and P-HM2)

### Contact Analysis ###
These scripts analyze protein interface contacts across the dark complex for CP12 interactions with other dark complex subunits GAPDH and PRK - and across variants and redox conditions. 
* **Data** is located in two directories ```Contact_Analysis/Binary/data``` (binary analysis, i.e. 0 if not or 1 if within cutoff, based on heavy atom cutoff) and ```Contact_Analysis/Contact_dynamics_analysis_current_withChemistryGrouping/data``` (comparative chemistry based analysis), these data are calculated from the trajectories available in trajectories.zip. 

#### Script Description and Guide ####
To perform contact analysis you first need to run analysis, then make a CSV from the output intermediate CSVs, and then subsequent figures such as heatmaps and other plots can be computed. This information is also available in ```/Contact_Analysis/README.md``` and each script header in the form of comments, but is more succintly outlined here:

STEP 1 - 
```analyze_cp12_contact_pairs.py``` - Computes residue-pair contact frequencies and differential contact rewiring maps comparing conditions. Takes as input molecular dynamics trajectories, structure files, and topologies and returns CSV files that record per-replicate pair frequencies, aggregated per-condition, collapsed to alignment coordinates, differential tables, top rewired pairs.

STEP 2 - ```interface_contacts_freq_create_master.py``` - Aggregates contact frequency data and generates a master CSV for visualization. 

#### Notes ####
The output from these above analysis and aggregation steps is used as input for a few different visualizations (see: ```./Contact_Analysis/visualizations```)

```./Contact_Analysis/Binary``` contains an individuation of each metric carried out by ```analyze_cp12_contact_pairs.py``` in aggregate

### MMGBSA ###
This directory contains the pipeline used to select, run, and extract Generalized-Born MM/GBSA
binding-energy calculations for the dark complex. 

It is necessary to have Gromacs installed in order to calculate the relative binding free energy. See `/MMGBSA/README.md` for more details. 


### SASA and BSA ###
This script calculates the solvent accessible surface area and the buried surface area, for different complex subunits. Note: a working GROMACS installation is required.  
* **Data** is located in ```SASA-BSA/data``` and contains a single .tsv recording the SASA and BSA for each subunit partner pair and the overall complex

### PCA ###
This subfolder contains the script along with the order for carrying out structural PCA on the dark complex. 

### MD_simulations ###
This subfolder contains the scripts needed to run molecular dynamics for all dark complex systems for 100ns in gromacs, and distinguishes between the 0000 (reduced) and 0003 (oxidized) simulations and the remaining 1010 PTM-modified systems. 

### PTM-plots ###

For more clear delineation each SI figure linked analysis is placed in a separate folder. 

```11``` - PCA distance plotted across system PTM number. 

```12``` - Log2 fold enrichment across each PTM type, for MED4 and P-HM2, and PTM ratio for each PTM type-specific site, colored by -log10 Benjamini-Hochberg

```13``` - Mean first passage time and residence time for the most conformationally distant systems. 

```14``` - Change in binding free energy as compared with each variants reference system (0003) for MED4 and P-HM2, all top 50 systems, and the 12 overlapping out of those. 

### RMSF-bin ###

Computes and plots CP12 Cα RMSF broken down by interface region (PRK-binding, GAP2-binding,
Non-interface) across all four dark-complex conditions (MED4 ox/red, HM2 ox/red) and three
replicates each.

### agent_skills_for_data_exploration ###

This contains agentic skills, prompt, and examples to help domain experts quickly explore MED4 proteomics data and triage for downstream analysis.
It requires VScode + Copilot and an LLM provider for the agent. See individual skills for their respective licenses.
* **Data** The processed MED4 proteomics files are available in the ```/data``` subfolder.* 

### Trajectory Data ###
In the parent Zenodo repository please find the ```trajectories.zip``` file containing a reference .pdb and three replicate .xtc files. These are broken into directories for MED4 and PHM2 systems, which contain a subfolder for each variant's 0000 system (fully reduced system with no modifications) and 0003 system (reference oxidized with disulfide modifications on CP12 and PRK subunits).

## Usage ##
Please refer to each sub-directory (listed above) for a full description of how to use the included scripts and accompanying data. 

## Provenance and Portability Note ##
This release is intended to preserve the original analysis provenance (scripts, parameters, and execution context) used in the study.
Several workflows retain environment-specific assumptions (for example HPC scheduler settings, filesystem layouts, and tool paths) and may require local adaptation before rerunning outside the original infrastructure.


## Authors ##
- Malio Nelson (First Author, Maintainer)
- Scott Widmann (First Author)
- Ruonan Wu (First Author)
- August George
- Xiaolu Li 
- Song Feng
- Margaret S. Cheung (Corresponding Author)
- Amity Andersen (Corresponding Author)
- Kehinde Idowu (Corresponding Author)


## Acknowledgements ##
This work was supported by the NW-BRaVE for Biopreparedness project funded by the U.S. Department of Energy (DOE), Office of Science, Office of Biological and Environmental Research, under FWP 81832. Pacific Northwest National Laboratory is a multi-program national laboratory operated by Battelle for the DOE under Contract DE-AC05-76RL01830. Portions of this research were performed at the Environmental Molecular Sciences Laboratory, a DOE Office of Science User Facility sponsored by the Biological and Environmental Research program, under Proposal 61054 (Award DOI: 10.46936/staf.proj.2023.61054/60012367). NW-BRaVE additionally leveraged resources from the Advanced Computing Ecosystem IRI testbed at Oak Ridge National Laboratory’s Leadership Computing Facility, allocations from the 2024 ASCR Leadership Computing Challenge (ALCC), and 2025 NASEP allocations for Doudna. The research described herein was funded by the Generative AI for Science, Energy, and Security Science & Technology Investment under the Laboratory Directed Research and Development Program at Pacific Northwest National Laboratory (PNNL). This work was also supported by the Center for AI and Center for Cloud Computing at PNNL. 



