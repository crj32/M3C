# M3C: Monte Carlo Reference-based Consensus Clustering

Genome-wide data is used to stratify large complex datasets into classes using class discovery algorithms. A widely applied technique is consensus clustering, however; the approach is prone to overfitting and false positives. These issues arise from not considering reference distributions while selecting the number of classes (K). As a solution, we developed Monte Carlo reference-based consensus clustering (M3C). M3C uses a multi-core enabled Monte Carlo simulation to generate null distributions along the range of K which are used to select its value. Using a reference, that maintains the correlation structure of the input features, eliminates the limitations of consensus clustering. M3C uses the Relative Cluster Stability Index (RCSI) and p values to decide on the value of K and reject the null hypothesis, K=1. M3C can also quantify structural relationships between clusters, and uses spectral clustering to deal with non-Gaussian and complex structures. M3C can automatically analyse biological or clinical data with respect to the discovered classes.  
  
Details:  
  
-M3C calculates the consensus rate, a measure of stability of the co-clustering of samples, which is quantified for each K using the PAC score  
-Generation of reference PAC distribution using a multi-core Monte Carlo simulation  
-Reference generation preserves feature-feature correlation structure of data  
-Using the reference distributions the Relative Cluster Stability Index (RCSI) and empirical p values are used to select K and reject the null, K=1.   
-Extrapolated p values are calculate by fitting a beta distribution  
-M3C with the RCSI has superior performance than other methods, it eliminates overfitting  
-Ability to investigates structural relationships using hierarchical clustering of medoids and Sigclust  
-Includes spectral clustering to investigate complex structures (e.g. non linear, anisotropic)  
-Automatic re ordering of expression matrix and annotation data to help user do their analysis faster  
-Automatic analysis of clinical or biological data using survival analysis, chi-squared or Kruskal-Wallis tests  
-Plotting code using ggplot2 for publication quality outputs  
-User friendly PCA and tSNE functions that interface with the results  

