# M3C: Monte Carlo Consensus Clustering

New consensus clustering algorithm that takes into account expectation under the null hypothesis using a monte carlo simulation, this is based on the original algorithm by Monti et al. (2003) and subsequent work by 
Șenbabaoğlu et al. (2014).

Algorithm:  
  
-Calculates the consensus rate, a measure of stability, which is quantified using PAC score  
-Generation of reference PAC distribution using monte carlo simulation  
-Randomisation preserves gene-gene correlation structure of data  
-The PAC statistic and empirical p values instead of delta K  
-Increased accuracy compared with other methods verified using simulations  
-Multi core enabled monte carlo  
-Controls for the null hypothesis K = 1  
-Extrapolated p values for estimation of p values below the min simulated PAC  
  
Other features:    
  
-Re ordering of expression matrix and annotation data to help user easily get results  
-ggplot2 plotting code for publicaton quality outputs    

References:  
  
Monti, Stefano, et al. "Consensus clustering: a resampling-based method for class discovery and visualization of gene expression microarray data." Machine learning 52.1-2 (2003): 91-118.

Șenbabaoğlu, Yasin, George Michailidis, and Jun Z. Li. "Critical limitations of consensus clustering in class discovery." Scientific reports 4 (2014).
