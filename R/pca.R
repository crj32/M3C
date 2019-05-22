#' pca: A principal component analysis function
#' 
#' This is a flexible PCA function that can be run on the results object from M3C or a standard data frame.
#' It is a wrapper for prcomp/ggplot2 code, and can be customised with different colours and font sizes.
#' 
#' @param mydata Data frame or matrix or M3C results object: if dataframe/matrix should have samples as columns and rows as features
#' @param printres Logical flag: whether to print the PCA into current directory
#' @param K Numerical value: if running on the M3C results object, which value was the optimal K?
#' @param labels Character vector: if we want to just label with gender for example
#' @param text Character vector: if we wanted to label the samples with text IDs to look for outliers
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' @param textlabelsize Numerical value: text inside plot label size 
#' @param legendtitle Character vector: text legend title   
#' @param controlscale Logical flag: whether to control the colour scale
#' @param scale Numerical value: 1=spectral palette, 2=manual low and high palette, 3=categorical labels
#' @param low Character vector: continuous scale low colour
#' @param high Character vector: continuous scale high colour
#' @param colvec Character vector: a series of colours in vector for categorical labels, e.g. c("sky blue", "gold")
#'
#' @return A PCA plot object
#' @export
#'
#' @examples
#' PCA <- pca(mydata)

pca <- function(mydata, K = FALSE, printres = FALSE, labels = FALSE, text = FALSE, axistextsize = 18,
                legendtextsize = 18, dotsize = 5, textlabelsize = 4, legendtitle = 'Group',
                controlscale = FALSE, scale = 1, low = 'blue', high = 'red', colvec = c("sky blue", "gold")){
  if (K == FALSE && labels == FALSE && text == FALSE){
    pca1 = prcomp(t(mydata))
    scores <- data.frame(pca1$x) # PC score matrix
    p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = dotsize) + 
      theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
      scale_colour_manual(values = c("1"='sky blue'))
    if (printres == TRUE){
      message('printing PCA to current directory...')
      png('PCApriorclustering.png', height = 20, width = 22, units = 'cm',
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
    pca1 = prcomp(t(mydata))
    scores <- data.frame(pca1$x) # PC score matrix
    p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = factor(annon$consensuscluster)), size = dotsize) + 
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
      message('printing PCA to current directory...')
      png('PCApostclustering.png', height = 20, width = 26, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }else if (K == FALSE && labels != FALSE){ ################ KEY
    pca1 = prcomp(t(mydata))
    scores <- data.frame(pca1$x) # PC score matrix
    if (controlscale == TRUE){
      if (scale == 1){
        p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")
        #scale_colour_gradient(low="red", high="white")
      }else if (scale == 2){
        p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + #scale_colour_distiller(palette = "Spectral")
        scale_colour_gradient(low=low, high=high)
      }else if (scale == 3){
        p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) +
          scale_colour_manual(values = colvec)
      }
    }else{
      p <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = labels), size = dotsize) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = axistextsize, colour = 'black'),
              axis.text.x = element_text(size = axistextsize, colour = 'black'),
              axis.title.x = element_text(size = axistextsize),
              axis.title.y = element_text(size = axistextsize),
              legend.title = element_text(size = legendtextsize),
              legend.text = element_text(size = legendtextsize)) + 
        #guides(colour=guide_legend(title=legendtitle)) +
        labs(colour = legendtitle)
    }
    if (printres == TRUE){
      message('printing PCA to current directory...')
      png('PCAlabeled.png', height = 20, width = 26, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }else if (K == FALSE && labels == FALSE && text != FALSE){
    pca1 = prcomp(t(mydata))
    scores <- data.frame(pca1$x) # PC score matrix
    scores$label <- text
    p <- ggplot(data = scores, aes(x = PC1, y = PC2, label = label) ) + 
      geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = dotsize) + 
      theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
      scale_colour_manual(values = c("1"='sky blue')) + geom_text(vjust="inward",hjust="inward",size=textlabelsize)
    if (printres == TRUE){
      message('printing PCA to current directory...')
      png('PCApriorclustering.png', height = 20, width = 22, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }else{
    message('no valid options detected')
  }
  return(p)
}