#' pca: A principal component analysis function
#'
#' @param mydata Data frame or matrix or M3C results object: if dataframe/matrix should have samples as columns and rows as features
#' @param printres Logical flag: whether to print the PCA into current directory
#' @param K Numerical value: if running on the M3C results object, which value was the optimal K?
#'
#' @return A PCA plot object
#' @export
#'
#' @examples
#' PCA <- pca(mydata)

pca <- function(mydata, K = FALSE, printres = FALSE){
  if (K == FALSE){
    pca1 = prcomp(t(mydata))
    scores <- data.frame(pca1$x) # PC score matrix
    p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = 6) + 
      theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 30, colour = 'black'),
            axis.text.x = element_text(size = 30, colour = 'black'),
            axis.title.x = element_text(size = 30),
            axis.title.y = element_text(size = 30)) +
      scale_colour_manual(values = c("1"='sky blue'))
    if (printres == TRUE){
      message('printing PCA to current directory...')
      png('PCApriorclustering.png', height = 20, width = 22, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }else{
    res <- mydata
    mydata <- res$realdataresults[[K]]$ordered_data
    annon <- res$realdataresults[[K]]$ordered_annotation
    annon$id <- row.names(annon)
    annon <- annon[match(colnames(mydata), annon$id),]
    pca1 = prcomp(t(mydata))
    scores <- data.frame(pca1$x) # PC score matrix
    p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = factor(annon$consensuscluster)), size = 6) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 30, colour = 'black'),
            axis.text.x = element_text(size = 30, colour = 'black'),
            axis.title.x = element_text(size = 30),
            axis.title.y = element_text(size = 30),
            legend.title = element_text(size = 30),
            legend.text = element_text(size = 20)) + 
      guides(colour=guide_legend(title="Cluster"))
    if (printres == TRUE){
      message('printing PCA to current directory...')
      png('PCApostclustering.png', height = 20, width = 26, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }
  return(p)
}
