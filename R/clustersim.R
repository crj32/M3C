#' clustersim: A cluster simulator for testing clustering algorithms
#'
#' @param n Numerical value: The number of samples, it must be square rootable
#' @param n2 Numerical value: The number of features
#' @param r Numerical value: The radius to define the initial circle (use approx n/100)
#' @param K Numerical value: How many clusters to simulate
#' @param alpha Numerical value: How far to pull apart the clusters
#' @param wobble Numerical value: The degree of noise to add to the sample co ordinates
#' @param print Logical flag: whether to print the PCA into current directory
#' @param seed Numerical value: fixes the seed if you want to repeat results
#'
#' @return A list: containing 1) matrix with simulated data in it
#' @export
#'
#' @examples
#' res <- clustersim(225, 900, 8, 4, 0.75, 0.025, print = TRUE, seed=123)

clustersim <- function(n, n2, r, K, alpha, wobble, print = FALSE, seed=NULL){
  
  message('***clustersim***')
  
  if (n %% sqrt(n) != 0){
    stop("n or number of samples must have an integer as the square root")
  }
  if (is.null(seed) == FALSE){
    set.seed(seed)
  }
  addnoise <- function(m) { # wobble samples in circle a little
    rnorm(1, m, wobble*m) # usually 0.025
  }
  # draw 2, n2 element vectors from (0,1) normal distribution
  x <- rnorm(n2, mean = 0, sd = 0.1)
  y <- rnorm(n2, mean = 0, sd = 0.1)
  
  # create the grid - usual distribution
  matrix = matrix(nrow = n, ncol = 3)
  matrix[,1] <- rep(c(1:sqrt(n)),each=sqrt(n))
  matrix[,2] <- rep(c(1:sqrt(n)), sqrt(n))
  
  # remove points outside of r size
  x1 <- (cluster::pam(data.frame(matrix)[,1:2], 1)$medoids)[1]
  y1 <- (cluster::pam(data.frame(matrix)[,1:2], 1)$medoids)[2]
  for (q in seq(1,n,1)){ # compute distance from origin
    x2 <- matrix[q,1]
    y2 <- matrix[q,2]
    answer <- sqrt(((x2-x1)*(x2-x1)) + ((y2-y1)*(y2-y1)))
    matrix[q,3] <- answer
  }
  
  matrix2 <- subset(data.frame(matrix), X3 < r)
  #plot(matrix2[,1], matrix2[,2])
  matrix2[] <- vapply(as.matrix(matrix2), addnoise, numeric(1))
  #plot(matrix2[,1], matrix2[,2])
  
  # start with eigenvectors (rnorm created), and co ordinates which represent PCs
  # and this loop will simulate the data (n,n2) dimensions
  message('computing simulated data using linear combination...')
  res = matrix(nrow = nrow(matrix2), ncol = n2) # create the simulated data in circle shape
  for (i in seq(1,nrow(matrix2),1)){
    a <- matrix2[i,1]
    b <- matrix2[i,2]
    for (k in seq(1,n2,1)){
      xk <- x[k]
      yk <- y[k]
      answer <- a*xk + b*yk
      res[i,k] <- answer
    }
  }

  # take res (your simulated circle) and pull apart the data using K clusters
  message('calculating PCs of data and performing kmeans...')
  mydata <- as.data.frame(res)
  pca1 = prcomp(mydata)
  scores <- data.frame(pca1$x) # PC score matrix
  kmeanscluster <- kmeans(scores, K, nstart = 20) # calculate k means of the PC co ordinates
  clusters_km <- kmeanscluster$cluster
  clusters_km_df <- data.frame(clusters_km)
  clusters_km_df$ID <- row.names(clusters_km_df)
  
  # now need to seperate out the clusters
  scores$ID <- row.names(scores)
  merged <- merge(scores, clusters_km_df, by = 'ID')
  merged2 <- merged[c('PC1', 'PC2', 'clusters_km')]
  merged2 <- merged2[with(merged2, order(clusters_km)), ]
  merged2$PC1shifted <- NA
  merged2$PC2shifted <- NA
  cluster_sizes <- as.data.frame(table(merged2$clusters_km))
  mylist <- list()
  message('computing centroids of clusters and pulling apart...')
  for (i in seq(1,K,1)){
    # get co ordinates of centroid of cluster i
    xxx <- subset(merged2, clusters_km == i)
    x <- (cluster::pam(xxx[,1:2], 1)$medoids)[1] # centroid x co ordinate
    y <- (cluster::pam(xxx[,1:2], 1)$medoids)[2] # centroid y co ordinate
    # manipulate data frame, xxx, for cluster i
    xxx$PC1shifted <- NA
    xxx$PC2shifted <- NA
    for (j in seq(1:nrow(xxx))){
      xxx[j,4] <- xxx[j,1]+(alpha*x) # shift x co ordinate
      xxx[j,5] <- xxx[j,2]+(alpha*y) # shift y co ordinate
    }
    mylist[[i]] <- xxx # add the transformed data for cluster i into mylist
    if (i == K){
      final_df <- do.call("rbind", mylist)
    }
  }
  
  # convert back to a matrix of data for clustering
  final_df$ID <- row.names(final_df)
  final_df$ID <- as.numeric(final_df$ID)
  final_df <- final_df[with(final_df, order(ID)), ]
  final_matrix <- as.matrix(final_df)
  
  # transform back to n dimensional dataset
  message('transforming pulled apart PC co ordinates back to dataset...')
  jjj <- t(final_matrix[,4:5] %*% t(pca1$rotation[,1:2])) + pca1$center # equation, PCs * eigenvectors = original data
  # take jjj and do a PCA just to check the new dataset
  mydata2 <- as.data.frame(jjj)
  pca1 = prcomp(t(mydata2))
  scores <- data.frame(pca1$x) # PC score matrix
  p <- ggplot2::ggplot(data = scores, ggplot2::aes(x = PC1, y = PC2) ) + ggplot2::geom_point(ggplot2::aes(colour = factor(final_df$clusters_km)), size = 6) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.position = "none", panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = 33, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = 33, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = 33),
          axis.title.y = ggplot2::element_text(size = 33)) # 31 old
  if (print == TRUE){ # width 22 for figure1, width xx for figure 2
    png(paste('PCAsim',n,n2,r,K,alpha,wobble,'.png'), height = 14, width = 22, units = 'cm', # 16 old
         res = 900, type="cairo")
  }
  print(p) # print ggplot CDF in main plotting window
  if (print == TRUE){
    dev.off()
  }
  if (print == TRUE){
    print(p)
  }
  message('finished.')
  outputs <- list(simulated_data = jjj)
  return(outputs)
}





