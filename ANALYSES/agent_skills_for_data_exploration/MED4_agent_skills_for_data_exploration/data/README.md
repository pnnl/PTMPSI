# Readme for Processed Proteomics

## Experiment

Diel-adapted MED4 multi-PTM proteomics data taken at dawn (from 6:00 am to 6:30 am, with light perturbation and controls) and dust (from 6:00 pm to 6:30 pm, with dark perturbation and controls). Experiment begins at 6AM in darkness. Samples are exposed to light beginning immediately after the first measurement, with treatment in light and control continue in dark until 6:30AM, and the light is removed immediately after the measurement at 6PM with treatment in dark and control continue in light until 6:30PM.

## Files

* `BRaVE_MED4_light-dark_dawn-dusk_all_ptms.tsv`
    * Normalized intensities and statistical information for all PTMs (redox and acetylation).
* `BRaVE_MED4_light-dark_dawn-dusk_all_ptms_not_abundance_corrected.tsv`
    * Normalized intensities and statistical information for all PTMs (redox and acetylation), but without corresponding protein abundance correction.
* `BRaVE_MED4_light-dark_dawn-dusk_protein_abundance.tsv`
    * Normalized intensities and statistical information for protein abundance.

## Column Information
<!--#### `BRaVE_MED4_light-dark_dawn-dusk_all_ptms.tsv`-->
* `Peptide` : Peptides that are associated with the PTM site
* `Site` : PTM Site name, normally is the pattern of [protein name]@[Residue name] the protein name is normally some unique protein ID, such as Uniprot ID or annotated ID from fasta file used for peptide searching. For instance G6PDH@C35 means protein G6PDH's cysteine residue with AA sequence number 35.
* `UniProt` : This is normally uniprot ID for some model organism or well-annotated organism who has the proteome in the UniProt database. If it is uncommon organism, this column are just whatever ID was provided in the fasta file of the proteome sequences.
* `Protein` : similarly, this is protein name from UniProt if uniprot proteome fasta file was provided, if not from UniProt the name is likely to be the same as the UniProt column's IDs in fasta file.
* `Residue` : the residue name of the PTM sites, some times when combined with global proteomics data, the abundance of protein are annotated with Abund.
* `id` : this is the unique id for the PTM site, usually it will be the same as the `Site`.
* `Total Missingness` : this is the number of missing value in the proteomics data, specifically in all the intensities columns, see following columns.
* `Dark_0600 missingness` to `Light_1830 missingness` are the count of missing values in all the conditions, each condition is [Light or Dark]_[Time in 24hr format].
* `missing_check` : this is to check how many conditions has more than half of the replicate samples are missing. This is normally used for filtering out the data, because when there are too many missing values, it is not possible or less confident to perform statistical comparisons.
* `Dark_0600_R1` to `Light_1830_R5` : these are the spectral intensity values for each sample, the sample name in the pattern of [Light or Dark conditions]_[Time point with 24HR format]_[Replicate number]. 
* `ANOVA_[XXXX]_pval` : these are ANOVA p values of each factors and their interactions. 
* `Light_0630/Dark_0600` : the fold change between two conditions `Light_0630` and `Dark_0600`
* `Light_0630/Dark_0600_pval` : the p value of the t-test between two conditions `Light_0630` and `Dark_0600`
* `Light_0630/Dark_0600_adj-p` : the adjusted p value of t-test between two conditions `Light_0630` and `Dark_0600`
* Similarly for the rest of the comparisons between the conditons
* `Type` : this column annotates the types of PTM, Ox means oxidation (Cysteine), Ac means acetylation (Lysine)
* `site_number` : is calcuating how many PTM sites are captured in the same protein.
