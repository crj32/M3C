#' featurefilter: A function for filtering features
#'
#' This function is to filter features based on the coefficient of variation (A) or its second
#' order derivative (A2) (Kvalseth, 2017). A is standardised according to the mean so is better
#' for gene/ feature expression than variance or standard deviation, however, it comes with 
#' disadvantages. This is why we have included its second order derivative.
#'
#' @param mydata Data frame: should have samples as columns and rows as features
#' @param percentile Numerical value: the top X percent most variable features should be kept
#' @param method Character vector: coefficient of variation (A), the A second order derivative (A2)
#' @param topN Numerical value: the number of most variable features to display
#'
#' @return A list, containing: 
#' 1) filtered data
#' 2) statistics for each feature order according to A or A2
#' 
#' @references 
#' Reference for second order coefficient of variation
#' Kvålseth, Tarald O. "Coefficient of variation: the second-order alternative." Journal of Applied Statistics 44.3 (2017): 402-415.
#' 
#' @export
#'
#' @examples
#' filtered <- featurefilter(mydata,percentile=10)

featurefilter <- function(mydata,percentile=10,method='A',topN=20){
  
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
    # calc coefficient of variation for all rows (features)
    CV <- sigma/u
    A <- CV
    # get features with CV in the given percentile
    CVthresh <- quantile(CV, percentile) 
  }else if (method == 'A2'){
    message('performing calculations for second order co efficient of variation/A2')
    # calculations
    u <- rowMeans(mydata)
    sigma <- apply(mydata,1,sd)
    A <- sigma/u
    AA <- A^2
    # get second order co efficient of variation
    A2 <- sqrt((AA/(AA+1)))
    CV <- A2
    # get features with CV in the given percentile
    CVthresh <- quantile(CV, percentile)
  }
  
  ## filter data
  names <- names(CV)[CV>=as.numeric(CVthresh)]
  filtered_data <- subset(mydata, row.names(mydata) %in% names)
  
  # make data frame of results
  if (method == 'A'){
    test <- data.frame('feature'=row.names(mydata),'mean'=u,'sd'=sigma,'A'=A)
    test <- test[order(-test[,4]), ]
    message('printing topN most variable features with statistics...')
    print(head(test,topN))
  }else if (method == 'A2'){
    test <- data.frame('feature'=row.names(mydata),'mean'=u,'sd'=sigma,'A'=A,'A2'=A2)
    test <- test[order(-test[,5]), ]
    message('printing topN most variable features with statistics...')
    print(head(test,topN))
  }
  
  ## make results list
  mylist <- list('filtered_data'=filtered_data,'statistics'=test)
  
  # message
  message(paste('features remaining:',nrow(filtered_data)))
  
  #
  return(mylist)
}