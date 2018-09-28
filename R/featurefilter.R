#' featurefilter: A function for filtering genes based on the coefficient of variance (CV)
#'
#' @param mydata Data frame or matrix or M3C results object: if dataframe/matrix should have samples as columns and rows as features
#' @param mypercentile Numerical value: CV percentile by which to filter the features, e.g. 90 indicates top 10% based on CV
#'
#' @return A filtered (fewer rows) data frame
#' @export
#'
#' @examples
#' filtered <- featurefilter(mydata,mypercentile=90)

featurefilter <- function(mydata,mypercentile=75){
  message(paste('extracting the most variable: '),mypercentile,'% of features')
  message(paste('features to start with:',nrow(mydata)))
  # mypercentile (convert to decimal below)
  mypercentile <- mypercentile/100
  # calc co efficient of variation for all rows (features)
  CV <- apply(mydata,1,sd)/(rowMeans(mydata))
  # get features with CV in the given percentile
  CVthresh <- quantile(CV, mypercentile) 
  # filter data
  names <- names(CV)[CV>=as.numeric(CVthresh)]
  filtered_data <- subset(mydata, row.names(mydata) %in% names)
  # message and output data
  message(paste('features remaining:',nrow(filtered_data)))
  return(filtered_data)
}



