#' umap: A umap function
#'
#' @param mydata Data frame or matrix or M3C results (list) object: if dataframe/matrix should have samples as columns and rows as features
#' @param printres Logical flag: whether to print the plot into current directory
#' @param K Numerical value: if running on the M3C results object, which value was the optimal K? Needs manual input from user.
#' @param labels Factor: if we want to just display gender for example, only for when running without K parameter and with a matrix or data frame
#' @param seed Numerical value: to repeat the results exactly, setting seed is required
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' @param textlabelsize Numerical value: text inside plot label size
#'
#' @return A umap plot object
#' @export
#'
#' @examples
#' UMAP <- umap(mydata)
umap <- function(mydata, K=FALSE, labels=FALSE, printres=FALSE, seed=FALSE, axistextsize = 18,
                 legendtextsize = 18, dotsize = 5, textlabelsize = 4){
  if (seed != FALSE){
    set.seed(seed)
  }
  if (K == FALSE && labels == FALSE){
    umap <- umap::umap(t(as.matrix(mydata)))
    scores <- data.frame(umap$layout) # PC score matrix
    p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = dotsize) +
      theme_bw() +
      theme(legend.position = "none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize))+
      scale_colour_manual(values = c("1"='sky blue'))
    if (printres == TRUE){
      message('printing UMAP to current directory...')
      png('UMAPpriorclustering.png', height = 20, width = 22, units = 'cm',
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
    umap <- umap::umap(t(as.matrix(mydata)))
    scores <- data.frame(umap$layout) # PC score matrix
    p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = factor(annon$consensuscluster)), size = dotsize) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize),
            legend.title = element_text(size = legendtextsize),
            legend.text = element_text(size = legendtextsize)) +
      guides(colour=guide_legend(title="Cluster"))
    if (printres == TRUE){
      message('printing UMAP to current directory...')
      png('UMAPpostclustering.png', height = 20, width = 26, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }else if (K == FALSE && labels != FALSE){
    umap <- umap::umap(t(as.matrix(mydata)))
    scores <- data.frame(umap$layout) # PC score matrix
    p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize),
            legend.title = element_text(size = legendtextsize),
            legend.text = element_text(size = legendtextsize)) +
      guides(colour=guide_legend(title="Cluster"))
    if (printres == TRUE){
      message('printing UMAP to current directory...')
      png('UMAPlabeled.png', height = 20, width = 26, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }
  return(p)
}