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
#' filtered <- featurefilter(mydata, percentile = 10)

featurefilter <- function(mydata, percentile = 10, method= 'MAD', topN = 20){
  
  message(paste0('Extracting the top ', percentile, '% most variable features'))
  
  message(paste('Features to start with:', nrow(mydata)))
  

# Convert Percent to Decimal ----------------------------------------------

  percentile <- 1 - (percentile/100)
  

# Calculate Mean and Variance ---------------------------------------------

  u <- rowMeans(mydata)
  sigma <- apply(mydata, 1, sd)
  vars <- sigma^2

  
# Perform Method-Specific Calculation -------------------------------------

  if (method == 'A'){
    message('Performing calculations for coefficient of variation (A)')
    # calc coefficient of variation for all rows (features)
    CV <- sigma/u
    CV[is.na(CV)] <- 0
  } else if (method == 'A2'){
    message('Performing calculations for second order coefficient of variation (A2)')
    CV <- sigma/u
    CV[is.na(CV)] <- 0
    A <- CV
    AA <- CV^2
    # get second order coefficient of variation
    A2 <- sqrt((AA/(AA + 1)))
    CV <- A2
  } else if (method == 'var'){
    message('Performing calculations for variance')
    CV <- vars
  } else if (method == 'MAD'){
    message('Performing calculations for median absolute deviation')
    MAD <- apply(mydata, 1, mad)
    CV <- MAD
  }
  

# Filter Data -------------------------------------------------------------

  # Get features with CV in the given percentile
  CVthresh <- quantile(CV, percentile, na.rm = TRUE)
  names <- names(CV)[CV >= as.numeric(CVthresh)]
  filtered_data <- subset(mydata, row.names(mydata) %in% names)
  

# Make Data Frame of Results ----------------------------------------------

  results <- data.frame('feature' = row.names(mydata), 'mean' = u, 'var' = vars, 'sd' = sigma)

  if (method == 'A'){
    results <- cbind(results, data.frame('A' = CV))
  } else if (method == 'A2'){
    results <- cbind(results, data.frame('A'= A, 'A2'= A2))
  } else if (method == 'var'){
  # Nothing needs to be done
  } else if (method == 'MAD'){
    results <- cbind(results, data.frame('MAD' = MAD))
  }
  
  results <- results[order(-results[, ncol(results)]), ]
  message('Printing topN most variable features with statistics...')
  print(head(results, topN))
  
  output <- list('filtered_data' = filtered_data, 'statistics' = results)

  message(paste('Features remaining:', nrow(filtered_data)))

  return(output)
}