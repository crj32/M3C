# M3C: Monte Carlo Reference-based Consensus Clustering

M3C is a consensus clustering algorithm that improves performance by eliminating overestimation of K and can test the null hypothesis K=1.  

Details:  
  
-M3C calculates the consensus rate, a measure of stability of the co-clustering of samples, for all samples for every K  
-Either the PAC score or entropy can be used to quantify how stable the consensus matrix is   
-Generation of null models using a multi-core Monte Carlo simulation    
-Reference generation preserves feature-feature correlation structure of data  
-Using the reference distributions the Relative Cluster Stability Index (RCSI) and empirical p values are used to select K and reject the null, K=1   
-Extrapolated p values are calculated by fitting normal or beta distributions  
-A second method is included for faster results that uses a penalty term, called regularised consensus clustering  
-Automatic re ordering of expression matrix and annotation data to help user do their analysis faster  
-User friendly PCA, tSNE, and UMAP functions    
-All plotting code using ggplot2 for publication quality outputs  
  
Usage:  
  
res <- M3C(mydata)   

