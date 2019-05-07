# M3C: Monte Carlo Reference-based Consensus Clustering

M3C is a consensus clustering algorithm that improves performance by eliminating overfitting, it can also test the null hypothesis K=1.  

Details:  
  
-M3C calculates the consensus rate, a measure of stability of the co-clustering of samples, which is quantified for each K using the PAC score  
-Generation of reference PAC distribution using a multi-core Monte Carlo simulation  
-Reference generation preserves feature-feature correlation structure of data  
-Using the reference distributions the Relative Cluster Stability Index (RCSI) and empirical p values are used to select K and reject the null, K=1.   
-Extrapolated p values are calculate by fitting a beta distribution  
-A second method is included for faster results that uses a penalty term instead of a Monte Carlo simulation to deal with overfitting  
-Automatic re ordering of expression matrix and annotation data to help user do their analysis faster  
-Automatic analysis of clinical or biological data using survival analysis, chi-squared, or Kruskal-Wallis tests  
-User friendly PCA, tSNE, and UMAP functions that interface with the results  
-All plotting code using ggplot2 for publication quality outputs  
  
Usage:  
  
res <- M3C(mydata)   

