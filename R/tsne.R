#' tsne: A tsne function
#'
#' @param mydata Data frame or matrix or M3C results (list) object: if dataframe/matrix should have samples as columns and rows as features
#' @param printres Logical flag: whether to print the plot into current directory
#' @param K Numerical value: if running on the M3C results object, which value was the optimal K? Needs manual input from user.
#' @param labels Factor: if we want to just display gender for example, only for when running without K parameter and with a matrix or data frame
#' @param seed Numerical value: to repeat the results exactly, setting seed is required
#' @param perplex Numerical value: this is the perplexity parameter for tsne, it usually requires adjusting for each dataset
#'
#' @return A tsne plot object
#' @export
#'
#' @examples
#' TSNE <- tsne(mydata,perplex=15)

tsne <- function(mydata, K=FALSE, labels=FALSE, perplex=15, printres=FALSE, seed=FALSE){
  if (seed != FALSE){
    set.seed(seed)
  }
  if (K == FALSE && labels == FALSE){
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=TRUE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix
    p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = 6) + 
      theme_bw() + 
      theme(legend.position = "none", 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 30, colour = 'black'),
            axis.text.x = element_text(size = 30, colour = 'black'),
            axis.title.x = element_text(size = 30),
            axis.title.y = element_text(size = 30))+
      scale_colour_manual(values = c("1"='sky blue'))
    if (printres == TRUE){
      message('printing tSNE to current directory...')
      png('TSNEpriorclustering.png', height = 20, width = 22, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }else if (K != FALSE && labels == FALSE){
    res <- mydata
    mydata <- res$realdataresults[[K]]$ordered_data
    annon <- res$realdataresults[[K]]$ordered_annotation
    annon$id <- row.names(annon)
    annon <- annon[match(colnames(mydata), annon$id),]
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=TRUE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix
    p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = factor(annon$consensuscluster)), size = 6) + 
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
      message('printing tSNE to current directory...')
      png('TSNEpostclustering.png', height = 20, width = 26, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }else if (K == FALSE && labels != FALSE){
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=TRUE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix
    p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = factor(labels)), size = 6) + 
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
      message('printing tSNE to current directory...')
      png('TSNElabeled.png', height = 20, width = 26, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }
  return(p)
}