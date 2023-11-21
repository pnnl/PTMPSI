## Clustering analysis

First use the notebooks Concatenate_PDB.ipynb and preprocessing_concatenated_PDB.ipynb to concatenate the target PDB trajectories and preprocess before analyzing. e.g. Discard the last column. 

Then run CalcD_permutation_max.py to generate the matrix data for clustering map plotting. An example command would be:

```
python CalcD_permutation_max.py example.pdb output.dat 2 1 0 4
```

Finally, with the output .dat data file, run seaborn_plot_cluster_map to plot the cluster map.
