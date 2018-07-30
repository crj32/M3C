# M3C: Monte Carlo Reference-based Consensus Clustering

Genome-wide data is used to stratify patients into classes using class discovery algorithms. A popular option is consensus clustering, however; the method is prone to overfitting and false positives. These arise from not considering reference distributions while selecting the number of classes (K). As a solution, we developed Monte Carlo reference-based consensus clustering (M3C). M3C uses a Monte Carlo simulation to generate null distributions along the range of K which are used to select its value. Including a reference removes the limitations of the consensus clustering method. M3C can also quantify structural relationships between clusters and uses self-tuning spectral clustering to deal with non-Gaussian and complex structures.

Details:  
  
-M3C calculates the consensus rate, a measure of stability of samples, which is quantified for each K using the PAC score  
-Generation of reference PAC distribution using a multi-core Monte Carlo simulation  
-Reference generation preserves feature-feature correlation structure of data  
-The Relative Cluster Stability Index (RCSI) and empirical p values are used to select K and reject the null, K=1.   
-Extrapolated p values are calculate by fitting a beta distribution  
-M3C with the RCSI has better accuracy than other methods, verified using simulations, and eliminates overfitting
-Ability to investigates structural relationships using hierarchical clustering of medoids and Sigclust  
-Includes spectral clustering to investigate complex structures (e.g. non linear, anisotropic)    
-Automatic re ordering of expression matrix and annotation data to help user do their analysis faster  
-Plotting code using ggplot2 for publication quality outputs  
-User friendly PCA and tSNE functions that interface with the results

