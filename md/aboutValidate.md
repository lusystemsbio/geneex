#### Validate

  The “Validate” tab is available for users to compare simulated expression data to some reference experimental or simulated expression data to investigate whether the GRC simulations capture the essential features of the experimental data. Simulated data can be uploaded by the user from a file or directly from the RACIPE tab following a simulation. 


The comparison protocol begins by clustering the reference data into a number of clusters specified by the user. Next, each simulated model (“simulated sample”) is compared to every reference sample in each cluster, and subsequently assigned to either one of the reference clusters/phenotypes or to a null cluster if its expression pattern is not consistent with any of the experimentally observed expression clusters. The null hypothesis is generated from a user-specified number of gene permutations, with the p-value cutoff being specified by the user as well. The percentage of simulated models and reference samples belonging to each cluster in the reference are reported alongside the overall Kullback–Leibler divergence between the two distributions. Heatmaps of the simulated and reference expression datasets, as well as of the sample-sample correlations, are also plotted. Cluster 0 is the null cluster and by default we add one sample belonging to null cluster to the reference samples.

<img src="images/validate_1.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:300px; height:300px" />

<br> 

