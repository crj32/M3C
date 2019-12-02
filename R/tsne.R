#' tsne: A t-SNE function
#'
#' This is a flexible t-SNE function that can be run on a standard data frame.
#' It is a wrapper for Rtsne/ggplot2 code and can be customised with different colours and font sizes and more.
#' 
#' @param mydata Data frame or matrix: if dataframe/matrix should have samples as columns and rows as features
#' @param printres Logical flag: whether to print the t-SNE into current directory
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
#' @param printheight Numerical value: png height
#' @param printwidth Numerical value: png width
#' @param seed Numerical value: optionally set the seed
#' @param perplex Numerical value: perplexity value that Rtsne uses internally
#'
#' @return A t-SNE plot object
#' @export
#'
#' @examples
#' TSNE <- tsne(mydata,perplex=15)

tsne <- function(mydata, labels=FALSE, perplex=15, printres=FALSE, seed=FALSE, axistextsize = 18,
                 legendtextsize = 18, dotsize = 5, textlabelsize = 4, legendtitle = 'Group',
                 controlscale = FALSE, scale = 1, low = 'grey', high = 'red', 
                 colvec = c("skyblue", "gold", "violet", "darkorchid", "slateblue", "forestgreen", 
                            "violetred", "orange", "midnightblue", "grey31", "black"),
                 printheight = 20, printwidth = 22, text = FALSE){
  
  ## basic error handling
  
  if ( controlscale == TRUE && class(labels) %in% c( "character", "factor") && scale %in% c(1,2) ) {
    stop("when categorical labels, use scale=3")
  }
  if ( controlscale == TRUE && class(labels) %in% c( "numeric") && scale %in% c(3) ) {
    stop("when continuous labels, use scale=1 or scale=2")
  }
  if ( controlscale == FALSE && scale %in% c(2,3) ) {
    warning("if your trying to control the scale, please set controlscale=TRUE")
  }
  if (sum(is.na(labels)) > 0 && class(labels) %in% c('character','factor')){
    warning("there is NA values in the labels vector, setting to unknown")
    labels <- as.character(labels)
    labels[is.na(labels)] <- 'Unknown'
  }
  if (sum(is.na(text)) > 0 && class(text) %in% c('character','factor')){
    warning("there is NA values in the text vector, setting to unknown")
    text <- as.character(text)
    text[is.na(text)] <- 'Unknown'
  }
  
  ##
  
  message('***t-SNE wrapper function***')
  message('running...')
  
  if (seed != FALSE){
    set.seed(seed)
  }
  
  ## combinations
  
  # K FALSE, labels FALSE, text FALSE
  # K TRUE, labels FALSE (when K TRUE just over ride everything)
  # K FALSE, labels TRUE, text ?
  # require the above, but with text or without so new entry
  # K FALSE, labels FALSE, text TRUE
  
  ##
  
  if (labels[1] == FALSE && text[1] == FALSE){
    
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix
    
    p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(colour='skyblue', size = dotsize) + 
      theme_bw() + 
      theme(legend.position = "none", 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize))+
      scale_colour_manual(values = colvec)
    
    if (printres == TRUE){
      message('printing t-SNE to current directory...')
      png('TSNE.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
    
  }else if (labels[1] != FALSE && text[1] == FALSE){ ##### KEY
    
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix
    
    if (controlscale == TRUE){
      if (scale == 1){
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) + 
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
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) + 
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
        p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) + 
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
      p <- ggplot(data = scores, aes(x = X1, y = X2) ) + geom_point(aes(colour = labels), size = dotsize) + 
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
      message('printing tSNE to current directory...')
      png('TSNElabeled.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
    
  }else if (labels[1] != FALSE && text[1] != FALSE){ ##### KEY
    
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y) # PC score matrix
    scores$label <- text
    
    if (controlscale == TRUE){
      if (scale == 1){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(size = axistextsize, colour = 'black'),
                axis.text.x = element_text(size = axistextsize, colour = 'black'),
                axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) + 
          #guides(colour=guide_legend(title=legendtitle)) +
          labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
        #scale_colour_gradient(low="red", high="white")
      }else if (scale == 2){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) + 
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
          scale_colour_gradient(low=low, high=high)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
      }else if (scale == 3){
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) + 
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
          scale_colour_manual(values = colvec)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
      }
    }else{
      p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + geom_point(aes(colour = labels), size = dotsize) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = axistextsize, colour = 'black'),
              axis.text.x = element_text(size = axistextsize, colour = 'black'),
              axis.title.x = element_text(size = axistextsize),
              axis.title.y = element_text(size = axistextsize),
              legend.title = element_text(size = legendtextsize),
              legend.text = element_text(size = legendtextsize)) + 
        #guides(colour=guide_legend(title=legendtitle)) +
        labs(colour = legendtitle)+ geom_text(vjust="inward",hjust="inward",size=textlabelsize)
    }
    
    if (printres == TRUE){
      message('printing t-SNE to current directory...')
      png('TSNElabeled.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
    
  }else if (labels[1] == FALSE && text[1] != FALSE){
    
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
    scores <- data.frame(tsne$Y)
    scores$label <- text
    
    p <- ggplot(data = scores, aes(x = X1, y = X2, label = label) ) + 
      geom_point(aes(colour = factor(rep(1, ncol(mydata)))), size = dotsize) + 
      theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = 'black'),
            axis.text.x = element_text(size = axistextsize, colour = 'black'),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
      scale_colour_manual(values = colvec) + geom_text(vjust="inward",hjust="inward",size=textlabelsize)
    
    if (printres == TRUE){
      message('printing t-SNE to current directory...')
      png('TSNE.png', height = printheight, width = printwidth, units = 'cm',
          res = 900, type = 'cairo')
      print(p) # print ggplot CDF in main plotting window
      dev.off()
    }
  }
  
  message('done.')
  
  return(p)
}