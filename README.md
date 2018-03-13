# M3C: Monte Carlo Consensus Clustering

Genome-wide data is used to stratify patients into classes using class discovery algorithms. However, we have observed systematic bias present in current state-of-the-art methods. This arises from not considering reference distributions while selecting the number of classes (K). As a solution, we developed a consensus clustering-based algorithm with a hypothesis testing framework called Monte Carlo consensus clustering (M3C). M3C uses a multi-core enabled Monte Carlo simulation to generate null distributions along the range of K which are used to calculate p values to select its value. P values beyond the limits of the simulation are estimated using a beta distribution. M3C can quantify structural relationships between clusters and uses spectral clustering to deal with non-gaussian and imbalanced structures.

Details:  
  
-M3C calculates the consensus rate, a measure of stability of samples, which is quantified for each K using the PAC score  
-Generation of reference PAC distribution using a multi-core Monte Carlo simulation  
-Reference generation preserves gene-gene correlation structure of data  
-The relative cluster stability index (RCSI) and empirical p values are used instead of delta K 
-Extrapolated p values are calculate by fitting a beta distribution
-Increased accuracy compared with other methods verified using simulations  
-Controls for the null hypothesis K = 1  
-Removes systematic bias
-Ability to investigates structural relationships using hierarchical clustering of medoids and sigclust
-Inner algorithms are PAM, K means, and spectral clustering
-Automatic re ordering of expression matrix and annotation data to help user do their analysis faster
-Plotting code using ggplot2 for publication quality outputs    

