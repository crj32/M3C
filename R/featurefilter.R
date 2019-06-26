#' featurefilter: A function for filtering features
#'
#' This function is to filter features based on variance. Depending on the data different
#' metrics will be more appropiate, simple variance is included if variance does not tend to
#' increase with the mean. There is also the median absolute deviation which is a more robust
#' metric than variance, this is preferable. The coefficient of variation (A) or its second
#' order derivative (A2) (Kvalseth, 2017) are also included which standardise the standard
#' deviation with respect to the mean. It is best to manually examine the mean-variance relationship of 
#' the data, for example, using the results from this function together with the qplot function 
#' from ggplot2.
#'
#' @param mydata Data frame: should have samples as columns and rows as features
#' @param percentile Numerical value: the top X percent most variable features should be kept
#' @param method Character vector: variance (var), coefficient of variation (A), second order A (A2), median absolute deviation (MAD)
#' @param topN Numerical value: the number of most variable features to display
#'
#' @return A list, containing: 
#' 1) filtered data
#' 2) statistics for each feature order according to the defined filtering metric
#' 
#' @references 
#' Kvålseth, Tarald O. "Coefficient of variation: the second-order alternative." Journal of Applied Statistics 44.3 (2017): 402-415.
#' 
#' @export
#'
#' @examples
#' filtered <- featurefilter(mydata,percentile=10)

featurefilter <- function(mydata,percentile=10,method='MAD',topN=20){
  
  message('***feature filter function***')
  
  message(paste('extracting the most variable: '),percentile,' percent')
  message(paste('features to start with:',nrow(mydata)))
  
  # percentile (convert to decimal below)
  percentile <- 1-(percentile/100)
  
  if (method == 'A'){
    message('performing calculations for co efficient of variation/A')
    # calculate mean and variance
    u <- rowMeans(mydata)
    sigma <- apply(mydata,1,sd)
    vars <- sigma^2
    # calc coefficient of variation for all rows (features)
    CV <- sigma/u
    CV[is.na(CV)] <- 0
    A <- CV
    # get features with CV in the given percentile
    CVthresh <- quantile(CV, percentile, na.rm = TRUE) 
  }else if (method == 'A2'){
    message('performing calculations for second order co efficient of variation/A2')
    # calculations
    u <- rowMeans(mydata)
    sigma <- apply(mydata,1,sd)
    vars <- sigma^2
    A <- sigma/u
    A[is.na(A)] <- 0
    AA <- A^2
    # get second order co efficient of variation
    A2 <- sqrt((AA/(AA+1)))
    CV <- A2
    # get features with CV in the given percentile
    CVthresh <- quantile(CV, percentile, na.rm = TRUE)
  }else if (method == 'var'){
    message('performing calculations for variance')
    u <- rowMeans(mydata)
    sigma <- apply(mydata,1,sd)
    vars <- sigma^2
    CV <- vars
    CVthresh <- quantile(CV, percentile, na.rm = TRUE)
  }else if (method == 'MAD'){
    message('performing calculations for median absolute deviation')
    u <- rowMeans(mydata)
    MAD <- apply(mydata,1,mad)
    sigma <- apply(mydata,1,sd)
    vars <- sigma^2
    CV <- MAD
    CVthresh <- quantile(CV, percentile, na.rm = TRUE)
  }
  
  ## filter data
  names <- names(CV)[CV>=as.numeric(CVthresh)]
  filtered_data <- subset(mydata, row.names(mydata) %in% names)
  
  # make data frame of results
  if (method == 'A'){
    test <- data.frame('feature'=row.names(mydata),'mean'=u,'var'=vars,'sd'=sigma,'A'=A)
    test <- test[order(-test[,5]), ]
    message('printing topN most variable features with statistics...')
    print(head(test,topN))
  }else if (method == 'A2'){
    test <- data.frame('feature'=row.names(mydata),'mean'=u,'var'=vars,'sd'=sigma,'A'=A,'A2'=A2)
    test <- test[order(-test[,6]), ]
    message('printing topN most variable features with statistics...')
    print(head(test,topN))
  }else if (method == 'var'){
    test <- data.frame('feature'=row.names(mydata),'mean'=u,'var'=vars,'sd'=sigma)
    test <- test[order(-test[,3]), ]
    message('printing topN most variable features with statistics...')
    print(head(test,topN))
  }else if (method == 'MAD'){
    test <- data.frame('feature'=row.names(mydata),'mean'=u,'var'=vars,'sd'=sigma,'MAD'=MAD)
    test <- test[order(-test[,4]), ]
    message('printing topN most variable features with statistics...')
    print(head(test,topN))
  }
  
  ## make results list
  mylist <- list('filtered_data'=filtered_data,'statistics'=test)
  
  ## message
  message(paste('features remaining:',nrow(filtered_data)))
  
  ##
  return(mylist)
}