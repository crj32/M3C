#' M3C: Monte Carlo Consensus Clustering
#'
#' @param mydata Data frame or matrix: Contains the data, with samples as columns and rows as features
#' @param montecarlo Logical flag: whether to run the monte carlo simulation or not (recommended: TRUE)
#' @param cores Numerical value: how many cores to split the monte carlo simulation over
#' @param iters Numerical value: how many monte carlo iterations to perform (default: 100, recommended: 100-1000)
#' @param maxK Numerical value: the maximum number of clusters to test for, K (default: 10)
#' @param des Data frame: contains annotation data for the input data for automatic reordering (optional)
#' @param ref_method Character string: refers to which reference method to use (recommended: leaving this alone)
#' @param repsref Numerical value: how many reps to use for the monte carlo reference data (suggest 100)
#' @param repsreal Numerical value: how many reps to use for the real data (recommended: 100)
#' @param clusteralg String: dictates which algorithm to use for M3C (recommended: leaving this alone)
#' @param distance String: dictates which distance metric to use for M3C (recommended: leaving this alone)
#' @param pacx1 Numerical value: The 1st x co-ordinate for calculating the pac score from the CDF (default: 0.1)
#' @param pacx2 Numerical value: The 2nd x co-ordinate for calculating the pac score from the CDF (default: 0.9)
#' @param printres Logical flag: whether to print all results into current directory
#' @param printheatmaps Logical flag: whether to print all the heatmaps into current directory
#' @param showheatmaps Logical flag: whether to show the heatmaps on screen (can be slow)
#' @param seed Numerical value: fixes the seed if you want to repeat results, set the seed to 123 for example here
#'
#' @return A list: containing 1) the stability results and 2) all the output data (in another list) 
#' (see vignette for more details on how to access)
#' @export
#'
#' @examples
#' res <- M3C(mydata, cores=1, iters=50, ref_method = 'reverse-pca', montecarlo = TRUE,printres = FALSE, 
#' maxK = 10, showheatmaps = FALSE, repsreal = 50, repsref = 50,printheatmaps = FALSE, seed = 123, des = desx)
M3C <- function(mydata, montecarlo = TRUE, cores = 16, iters = 100, maxK = 10,
                              des = NULL, ref_method = 'reverse-pca', repsref = 100, repsreal = 100,
                              clusteralg = 'pam', distance = 'euclidean', pacx1 = 0.1, pacx2 = 0.9, printres = FALSE,
                              printheatmaps = FALSE, showheatmaps = FALSE, seed=NULL){

  if (is.null(seed) == FALSE){
    set.seed(seed)
  }

  message('***M3C: Monte Carlo Consensus Clustering***')

  # error handling of input variables

  if ( ! class( mydata ) %in% c( "data.frame", "matrix", "ExpressionSet" ) ) {
    stop("mydata must be a data frame, matrix, or ExpressionSet (eset object)")
  }
  if ( inherits( mydata,"ExpressionSet" ) ) {
    mydata <- exprs(mydata)
  }
  if (montecarlo != TRUE){
    message('warning running without monte carlo simulation lowers accuracy')
  }
  if (clusteralg != 'pam'){
    message('warning pam is strongly advisable, because it is far faster than kmeans and far more accurate than hc')
  }
  if (clusteralg == 'pam' && distance != 'euclidean'){
    message('warning pam must be used with euclidean distance, changing to euclidean...')
    distance <- 'euclidean'
  }
  if (clusteralg == 'km' && distance != 'euclidean'){
    message('warning kmeans must be used with euclidean distance, changing to euclidean...')
    distance <- 'euclidean'
  }
  if (ncol(mydata) > nrow(mydata)){
    message('samples(columns) exceeds features(rows), switching to cholesky decomposition method...')
    ref_method = 'chol' # this works when variables < samples
  }
  if (is.null(des) == TRUE){ # user does not give extra annotation data
    message('user does not give extra annotation for reordering according to clustering results')
  }
  if (is.null(des) == FALSE){ # user does give extra annotation data
    message('user does give extra annotation for reordering according to clustering results')
    '%!in%' <- function(x,y)!('%in%'(x,y))
    if("ID" %!in% colnames(des))
    {
      stop('in the supplied annotation object for reordering, ther is no \'ID\' column....')
    }
    if (all(colnames(mydata) %in% des$ID) == FALSE){
      stop('the IDs in the annotation do not match column IDs in data')
    }
    if (nrow(des) != ncol(mydata)){
      stop('the dimensions of your annotation object do not match data object')
    }
  }

  # consensuscluster2 reference or no reference functions

  if (montecarlo == TRUE){
    ## run monte carlo simulation to generate references with same gene-gene correlation structure
    message('running simulations...')
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    invisible(capture.output(pb <- txtProgressBar(min = 0, max = iters, style = 3)))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    ## for each loop to use all cores
    ls<-foreach::foreach(i = 1:iters, .export=c("ccRun", "CDF", "connectivityMatrix", "ConsensusClusterRef", "myPal",
                                       "sampleCols", "setClusterColors", "triangle", "rmv", "msr"),
                .packages=c("cluster", "base", "Matrix"), .combine='rbind', .options.snow = opts) %dopar% {

                  if (is.null(seed) == FALSE){
                    set.seed(i)
                  }

                  if (ref_method == 'reverse-pca'){ # reverse PCA method
                    pca1 = prcomp(t(mydata))
                    simulated_data <- matrix(nrow = ncol(mydata), ncol = ncol(mydata)) # make up the principle components
                    for (i in seq(1,ncol(pca1$x))){
                      simulated_data[,i] <- rnorm(nrow(pca1$x), mean = 0, sd = sd((pca1$x)[,i])) # this is the PC co ordinates
                    }
                    null_dataset <- t(t(simulated_data %*% t(pca1$rotation)) + pca1$center) # go back to null distrib
                    mydata2 <- t(as.data.frame(null_dataset))
                  }

                  if (ref_method == 'chol'){ # cholesky method
                    sorted <- nearPD(cov(t(mydata))) # converting to positive definate cov matrix
                    sorted2 <- rmv(305, sorted$mat, rfunc = rnorm, method = c("chol"))
                    sorted3 <- as.data.frame(t(as.matrix(sorted2)))
                    mydata2 <- sorted3
                  }

                  if (ref_method == 'within-range'){ # generate random data within the range of each variable
                    random <- matrix(nrow = nrow(mydata), ncol = ncol(mydata))
                    temp <- as.matrix(mydata)
                    for (jj in seq(1:nrow(temp))){
                      minv <- min(temp[jj,])
                      maxv <- max(temp[jj,])
                      random[jj,] <- runif(length(temp[jj,]), minv, maxv)
                    }
                    mydata2 <- as.data.frame(random)
                  }

                  m_matrix <- as.matrix(mydata2)
                  results <- ConsensusClusterRef(m_matrix,maxK=maxK,reps=repsref,pItem=0.8,pFeature=1, # only use 100x iterations
                                                 clusterAlg=clusteralg, # use pam it is fast
                                                 distance=distance, # with pam always use euclidean
                                                 title = '/home/christopher/Desktop/',
                                                 x1=pacx1, x2=pacx2, printres=FALSE, seed=seed)
                  pacresults <- results$pac_scores$PAC_SCORE
                  return(pacresults)
                }
    close(pb)
    parallel::stopCluster(cl)
    message('finished generating reference distribution')
    ## finished monte carlo, results are in ls matrix
    ## real data PAC score calculation
    results2 <- ConsensusClusterReal(as.matrix(mydata),maxK=maxK,reps=repsreal,pItem=0.8,pFeature=1, # only use 100x iterations
                                     clusterAlg=clusteralg, # use pam it is fast
                                     distance=distance, # with pam always use euclidean
                                     title = '/home/christopher/Desktop/',
                                     printres = printres,
                                     showheatmaps = showheatmaps, printheatmaps = printheatmaps, des = des,
                                     x1=pacx1, x2=pacx2, seed=seed) # png to file
    real <- results2$pac_scores
    allresults <- results2$allresults
    copres <- results2$copres # alternative metric

    ## process reference data and calculate scores and p values (simulations are in ls matrix)
    colnames(real)[2] <- 'PAC_REAL'
    real$PAC_REF <- colMeans(ls)
    real$PAC_STAT <- real$PAC_REF - real$PAC_REAL # alternative pac statistic without taking log

    # need to log it because pac real can reach 0
    pacreal <- real$PAC_REAL
    pacreal[pacreal==0] <- 0.01 # set zero really small
    PACREALLOG <- log(pacreal)
    PACREFLOG <- log(real$PAC_REF)
    real$PAC_STAT <- PACREFLOG - PACREALLOG # replace non log with log

    ## usual p value derivation
    pvals <- c()
    x = 1
    for (var in real$PAC_REAL){
      distribution <- as.numeric(ls[,x])
      pvals <- c(pvals,((length(distribution[distribution < real$PAC_REAL[x]])) + 1)/ (iters+1)) # (b+1)/(m+1)=pval
      x = x + 1
    }
    real$P_VALUE <- pvals # this object contains all the results
    real$P_SCORE <- -log10(real$P_VALUE)

    ## extrapolated p values, only use when observed PAC values is lower than all simulated PAC values
    ## really needs 1000 simulations
    ## using 10 lowest PAC scores to fit model to
    ## https://genepi.qimr.edu.au/staff/davidD/Sib-pair/Documents/extrapolatedMCPvalues2011.pdf (references and maths in here)
    extraps <- c() # results
    indices <- c() # store indices of p values to replace with extrapolated
    x = 1
    m = 10 # between 10 and 20 (number of top values used to calc)
    for (var in real$PAC_REAL){
      if (var < min(as.numeric(ls[,x]))){
        indices <- c(indices, x)
        distribution <- as.numeric(ls[,x])
        # invert
        distribution = 1 - distribution
        #
        lowestpacs <- sort(distribution, decreasing = TRUE)[1:(m+1)]
        summedvals <- 0
        finalval <- lowestpacs[(m+1)]
        for (val in lowestpacs[1:m]){
          ans <- log(val) - log(finalval)
          summedvals = summedvals + ans
        }
        pacobs <- real$PAC_REAL[x]
        # invert
        pacobs = 1 - pacobs
        #
        astar = (1/m)*summedvals
        pextra = (m/iters)*(pacobs/finalval)^(-1/astar)
        extraps <- c(extraps, pextra)
      }
      if (var >= min(as.numeric(ls[,x]))){
        extraps <- c(extraps, NA)
      }
      x = x + 1
    }
    # don't replace actual p values with extrapolated, but code is here
    real$P_EXTRAPOLATED <- extraps # attach extrapolated p values to results
    #real$P_VALUE[indices] <- extraps[indices]
    #real$P_SCORE <- -log10(real$P_VALUE)

    # plot real vs reference results
    # pac statistic
    px <- ggplot2::ggplot(data=real, ggplot2::aes(x=K, y=PAC_STAT)) + ggplot2::geom_line(colour = "purple", size = 2) + ggplot2::geom_point(colour = "black", size = 3) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 26, colour = 'black'),
            axis.text.x = ggplot2::element_text(size = 26, colour = 'black'),
            axis.title.x = ggplot2::element_text(size = 26),
            axis.title.y = ggplot2::element_text(size = 26),
            legend.text = ggplot2::element_text(size = 26),
            legend.title = ggplot2::element_text(size = 26),
            plot.title = ggplot2::element_text(size = 26, colour = 'black', hjust = 0.5),
            panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_x_continuous(breaks=c(seq(0,maxK,1))) +
      ggplot2::ylab('PAC Statistic') +
      ggplot2::xlab('Number of Clusters') +
      ggplot2::labs(title = "PAC Ref - PAC Real")

    # pval score
    col = ifelse(real$P_SCORE > 1.30103,'tomato','black')
    py <- ggplot2::ggplot(data=real, ggplot2::aes(x=K, y=P_SCORE)) + ggplot2::geom_point(colour = col, size = 3) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 26, colour = 'black'),
            axis.text.x = ggplot2::element_text(size = 26, colour = 'black'),
            axis.title.x = ggplot2::element_text(size = 26),
            axis.title.y = ggplot2::element_text(size = 26),
            legend.text = ggplot2::element_text(size = 26),
            legend.title = ggplot2::element_text(size = 26),
            plot.title = ggplot2::element_text(size = 26, colour = 'black', hjust = 0.5),
            panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_x_continuous(breaks=c(seq(0,maxK,1))) +
      ggplot2::ylab(expression('-log'[10]*'p')) +
      ggplot2::xlab('Number of Clusters') +
      ggplot2::labs(title = "PAC p values") +
      ggplot2::geom_hline(yintercept=1.30103, size=0.75, linetype='dashed', colour='tomato') # 0.05 sig threshold

    if (printres == TRUE){
      png(paste('pscore.png'), height = 14, width = 20, units = 'cm',
           res = 900, type = 'cairo')
    }
    print(py) # print ggplot CDF in main plotting window
    if (printres == TRUE){
      dev.off()
    }
    if (printres == TRUE){
      png(paste('pacdiff.png'), height = 14, width = 20, units = 'cm',
          res = 900, type = 'cairo')
    }
    print(px) # print ggplot CDF in main plotting window
    if (printres == TRUE){
      dev.off()
    }
    if (printres == TRUE){
      print(py)
      print(px)
    }
  }

  if (montecarlo == FALSE){
    message('running without monte carlo simulations...')
    results2 <- ConsensusClusterReal(as.matrix(mydata),maxK=maxK,reps=repsreal,pItem=0.8,pFeature=1, # pfeature = 1 default
                                     clusterAlg=clusteralg, # default = pam, others = hc, km
                                     distance=distance, # seed=1262118388.71279,
                                     title = '/home/christopher/Desktop/',
                                     printres = printres, x1=pacx1, x2=pacx2,
                                     showheatmaps = showheatmaps, printheatmaps = printheatmaps, des = des,
                                     seed=seed) # png to file
    real <- results2$pac_scores
    allresults <- results2$allresults
    copres <- results2$copres # alternative metric
  }

  if (file.exists('Rplots.pdf') == TRUE){ # random file pops up
    file.remove('Rplots.pdf') # this needs to be removed
    unlink('Rplots.pdf') # this needs to be removed
  }

  if (printres == TRUE){
    write.csv(real, file = 'pacresultfile.csv', row.names = FALSE)
  }

  # calinski code was here (in version18)
  # close(cl)

  # make an automated decision if running the monte carlo distribution
  if (montecarlo == TRUE){
    real$decision <- 'sub optimal'
    if (length(which(signif(real$P_VALUE, digits = 7) == min(signif(real$P_VALUE, digits = 7)))) > 1){ # multiple top hits
      message('we have multiple lowest p values, using pac stat to break the tie')
      pacstats <- real$PAC_STAT[(which(signif(real$P_VALUE, digits = 7) == min(signif(real$P_VALUE, digits = 7))))]
      real$decision[which(real$PAC_STAT == pacstats[(which(signif(pacstats, digits = 7) == max(signif(pacstats, digits = 7))))])] <- 'optimal'
      if (length(  (which(signif(pacstats, digits = 7) == max(signif(pacstats, digits = 7))))) > 2){
        message('multiple lowest p values and multiple highest pac stats, consider adding more iterations')
        real$decision <- 'tie'
      }
    }
    if (length(which(signif(real$P_VALUE, digits = 7) == min(signif(real$P_VALUE, digits = 7)))) == 1){ # just one top hit
      real$decision[which(signif(real$P_VALUE, digits = 7) == min(signif(real$P_VALUE, digits = 7)))] <- 'optimal'
    }
    # print the decision
    if (length(which(signif(real$P_VALUE, digits = 7) == min(signif(real$P_VALUE, digits = 7)))) > 1){ # multiple top hits
      message(c('optimal decision: ', which(real$decision == 'optimal')+1))
      if (length(  (which(signif(pacstats, digits = 7) == max(signif(pacstats, digits = 7))))) > 2){
        message('optimal decision: none')
      }
    }
    if (length(which(signif(real$P_VALUE, digits = 7) == min(signif(real$P_VALUE, digits = 7)))) == 1){ # just one top hit
      message(c('optimal decision: ', which(real$decision == 'optimal')+1))
    }
  }

  # return results with or without monte carlo
  return(list("allresults" = allresults, 'pac_scores' = real)) # contains pac scores, and cc matrices, reordered data + annotation
}

ConsensusClusterReal <- function( d=NULL, # function for real data
                                  maxK = 3,
                                  reps=10,
                                  pItem=0.8,
                                  pFeature=1,
                                  clusterAlg="hc",
                                  title="untitled_consensus_cluster",
                                  innerLinkage="average",
                                  finalLinkage="average",
                                  distance="pearson",
                                  ml=NULL,
                                  tmyPal=NULL,
                                  seed=NULL,
                                  weightsItem=NULL,
                                  weightsFeature=NULL,
                                  verbose=FALSE,
                                  corUse="everything",
                                  showheatmaps=FALSE,
                                  printheatmaps=FALSE,
                                  printres=FALSE,
                                  x1=0.1,
                                  x2=0.9,
                                  des = NULL) {
    message('running consensus cluster algorithm for real data...') # this is the main function that takes the vast majority of the time

    if (is.null(seed) == FALSE){
      set.seed(seed)
    }

    ml <- ccRun( d=d,
                 maxK=maxK,
                 repCount=reps,
                 diss=inherits(d,"dist"),
                 pItem=pItem,
                 pFeature=pFeature,
                 innerLinkage=innerLinkage,
                 clusterAlg=clusterAlg,
                 weightsFeature=weightsFeature,
                 weightsItem=weightsItem,
                 distance=distance,
                 verbose=verbose,
                 corUse=corUse)
    message('finished')
  res=list();
  colorList=list()
  colorM = rbind() #matrix of colors.
  #18 colors for marking different clusters
  thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
               "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776", "#ffffff")
  colBreaks = NA
  if (is.null(tmyPal) == TRUE) {
    colBreaks = 10
    tmyPal = myPal(colBreaks)
  }
  else {
    colBreaks = length(tmyPal)
  }
  sc = cbind(seq(0, 1, by = 1/(colBreaks)))
  rownames(sc) = sc[, 1]
  sc = cbind(sc, sc)

  resultslist <- list() # holds all results for each value of K

  for (tk in 2:maxK){
    fm = ml[[tk]]
    hc=hclust( as.dist( 1 - fm ), method=finalLinkage);
    ct = cutree(hc,tk)
    names(ct) = colnames(d)
    if(class(d)=="dist"){
      names(ct) = colnames(as.matrix(d))
    }
    c = fm
    colorList = setClusterColors(res[[tk-1]][[3]],ct,thisPal,colorList)
    pc = c

    pc=pc[hc$order,] #pc is matrix for plotting, same as c but is row-ordered and has names and extra row of zeros.
    pc = rbind(pc,0)
    colcols <- as.factor(as.numeric(as.factor(colorList[[1]]))) # this is for aheatmap/ other non heatmap engines
    cols <- colorRampPalette(RColorBrewer::brewer.pal(9,'Reds')[1:6])(256)
    if (showheatmaps == TRUE & printheatmaps == FALSE){
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
    }
    if (showheatmaps == TRUE & printheatmaps == TRUE){
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      png(paste('K=',tk,'heatmap.png'), height = 12, width = 12, units = 'cm',
          res = 900, type = 'cairo')
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      dev.off()
    }
    if (showheatmaps == FALSE & printheatmaps == TRUE){
      png(paste('K=',tk,'heatmap.png'), height = 12, width = 12, units = 'cm',
          res = 900, type = 'cairo')
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      dev.off()
    }
    res[[tk]] = list(consensusMatrix=c,consensusTree=hc,consensusClass=ct,ml=ml[[tk]],clrs=colorList)
    colorM = rbind(colorM,colorList[[1]])

    # fixed consensus matrix, now in correct order/ same order as reordered data
    usethisorder <- hc$order # is this correct
    colnames(pc) <- seq(1,ncol(pc))
    pc <- pc[,usethisorder] # get the consensus matrices in the correct order

    # code here for extracting ordered data for user
    m_matrix=as.matrix(d) # this is the input data which we will reorder to match the consensus clusters
    #cc_matrix <- pc # pc is the consensus matrix
    #colnames(cc_matrix) <- colnames(m_matrix)
    # o <- NMF::aheatmap(cc_matrix, scale = 'none', distfun = 'pearson', Rowv = NA, color = tmyPal,
    #               legend = FALSE, annLegend = FALSE, Colv = as.dendrogram(hc))
    clustering <- as.numeric(as.factor(colorList[[1]]))
    clusteringdf <- data.frame(sample = colnames(m_matrix), cluster = clustering)
    neworder1 <- colnames(m_matrix)[hc$order] # changed o$colInd to hc$order ### changed here
    neworder2 <- gsub('-', '.', neworder1)
    df <- data.frame(m_matrix)
    newdes <- data.frame(ID = colnames(df), consensuscluster = factor(clusteringdf$cluster))
    colnames(df) <- gsub('X', '', colnames(df))
    neworder2 <- gsub('X', '', neworder2)
    data <- df[neworder2]
    if (is.null(des) == TRUE){
      neworder1 <- gsub('-', '.', neworder1) # cc code is changing - to . so change back
      vec <- grepl('X', colnames(d)) # check for X's in original colnames if dont exist run this code
      if (all(vec == FALSE)){ # just changed this to FALSE
        newdes$ID <- gsub('X', '', newdes$ID) # this creates problem if X in original names
      }
      newerdes <- newdes[match(neworder1, newdes$ID),]
      annotation <- data.frame(newerdes$consensuscluster)
      row.names(annotation) <- newerdes$ID
      colnames(annotation) <- c('CONSENSUS CLUSTER')
    }
    if (is.null(des) == FALSE){
      neworder1 <- gsub('-', '.', neworder1)
      neworder1 <- gsub('X', '', neworder1)
      des$ID <- gsub('-', '.', des$ID) # this is a problem with formatting of the ids - check out later
      vec <- grepl('X', colnames(d)) # check for X's in original colnames if dont exist run this code
      if (all(vec == FALSE)){ # changed to true** then back to false, problems here
        newdes$ID <- gsub('X', '', newdes$ID) # this creates problem if X in original names
      }
      merged <- merge(newdes, des, by = 'ID') # name to merge by
      newdes <- merged
      newerdes <- newdes[match(neworder1, newdes$ID),]
      annotation <- newerdes
      row.names(annotation) <- newerdes$ID
      annotation$ID <- NULL
    }
    newList <- list("consensus_matrix" = pc, 'ordered_data' = data, 'ordered_annotation' = annotation) # you can remove ml
    resultslist[[tk]] <- newList
  }
  pac_res <- CDF(ml, printres=printres, x1=x1, x2=x2) # this runs the new CDF function with PAC score
  copres <- COP(ml) # this runs the new COP function
  res[[1]] = colorM
  listxyx <- list("allresults" = resultslist, 'pac_scores' = pac_res, 'copres' = copres) # use a name list, one item is a list of results
  return(listxyx)

  print('finished this function')

}

ConsensusClusterRef <- function( d=NULL, # function for reference data
                                 maxK = 3,
                                 reps=10,
                                 pItem=0.8,
                                 pFeature=1,
                                 clusterAlg="hc",
                                 title="untitled_consensus_cluster",
                                 innerLinkage="average",
                                 finalLinkage="average",
                                 distance="pearson",
                                 ml=NULL,
                                 tmyPal=NULL,
                                 weightsItem=NULL,
                                 weightsFeature=NULL,
                                 corUse="everything",
                                 x1=0.1,
                                 x2=0.9,
                                 printres=FALSE,
                                 seed=NULL) {
  if (is.null(seed) == FALSE){
    set.seed(seed)
  }
  verbose=FALSE
  showheatmaps=FALSE
  printheatmaps=FALSE
  printres=FALSE
  message('running consensus cluster algorithm for reference data...') # this is the main function that takes the vast majority of the time
    ml <- ccRun( d=d,
                 maxK=maxK,
                 repCount=reps,
                 diss=inherits(d,"dist"),
                 pItem=pItem,
                 pFeature=pFeature,
                 innerLinkage=innerLinkage,
                 clusterAlg=clusterAlg,
                 weightsFeature=weightsFeature,
                 weightsItem=weightsItem,
                 distance=distance,
                 verbose=verbose,
                 corUse=corUse)
    message('finished.')
  res=list();
  colorList=list()
  colorM = rbind() #matrix of colors.
  #18 colors for marking different clusters
  thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
               "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776", "#ffffff")
  for (tk in 2:maxK){
    if(verbose){
      message(paste("consensus ",tk))
    }
    fm = ml[[tk]]
    hc=hclust( as.dist( 1 - fm ), method=finalLinkage);
    ct = cutree(hc,tk)
    names(ct) = colnames(d)
    if(class(d)=="dist"){
      names(ct) = colnames(as.matrix(d))
    }
    c = fm
    colorList = setClusterColors(res[[tk-1]][[3]],ct,thisPal,colorList)
    pc = c
    pc=pc[hc$order,] #pc is matrix for plotting, same as c but is row-ordered and has names and extra row of zeros.
    pc = rbind(pc,0)
    colcols <- as.factor(as.numeric(as.factor(colorList[[1]]))) # CHRIS: this is for aheatmap/ other non heatmap engines
    cols <- colorRampPalette(RColorBrewer::brewer.pal(9,'Reds')[1:6])(256)
    if (showheatmaps == TRUE & printheatmaps == FALSE){
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
    }
    if (showheatmaps == TRUE & printheatmaps == TRUE){
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      tiff(paste('K=',tk,'heatmap.tiff'), height = 12, width = 12, units = 'cm',
           compression = "lzw", res = 600)
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      dev.off()
    }
    if (showheatmaps == FALSE & printheatmaps == TRUE){
      tiff(paste('K=',tk,'heatmap.tiff'), height = 12, width = 12, units = 'cm',
           compression = "lzw", res = 600)
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      dev.off()
    }
    res[[tk]] = list(consensusMatrix=c,consensusTree=hc,consensusClass=ct,ml=ml[[tk]],clrs=colorList)
    colorM = rbind(colorM,colorList[[1]])
  }
  pac_res <- CDF(ml, printres=FALSE, x1=x1, x2=x2) # this runs the new CDF function with PAC score
  res[[1]] = colorM
  newList <- list("consensus_matrices" = res, 'pac_scores' = pac_res) # now returning a fairly coherant list of results
  return(newList)
}

ccRun <- function( d=d,
                   maxK=NULL,
                   repCount=NULL,
                   diss=inherits( d, "dist" ),
                   pItem=NULL,
                   pFeature=NULL,
                   innerLinkage=NULL,
                   distance=NULL, #ifelse( inherits(d,"dist"), attr( d, "method" ), "euclidean" ),@@@@@
                   clusterAlg=NULL,
                   weightsItem=NULL,
                   weightsFeature=NULL,
                   verbose=NULL,
                   corUse=NULL) {
  m = vector(mode='list', repCount)
  ml = vector(mode="list",maxK)
  n <- ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  mCount = mConsist = matrix(c(0),ncol=n,nrow=n)
  ml[[1]] = c(0);
  if (is.null( distance ) ) distance <- 'euclidean'  ## necessary if d is a dist object and attr( d, "method" ) == NULL
  acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                            "pearson", "spearman" )
  main.dist.obj <- NULL
  if ( diss ){
    main.dist.obj <- d
    ## reset the pFeature & weightsFeature params if they've been set (irrelevant if d is a dist matrix)
    if ( ( !is.null(pFeature) ) &&
         ( pFeature < 1 ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n" )
      pFeature <- 1 # set it to 1 to avoid problems with sampleCols
    }
    if ( ! is.null( weightsFeature ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n" )
      weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
    }
  } else { ## d is a data matrix
    ## we're not sampling over the features
    if ( ( clusterAlg != "km" ) &&
         ( is.null( pFeature ) ||
           ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) ) {
      ## only generate a main.dist.object IFF 1) d is a matrix, 2) we're not sampling the features, and 3) the algorithm isn't 'km'
      if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=TRUE))!="function") ) stop("unsupported distance.")

        if(distance=="pearson" | distance=="spearman"){
          main.dist.obj <- as.dist( 1-cor(d,method=distance,use=corUse ))
        }else if( class(try(get(distance),silent=TRUE))=="function"){
          main.dist.obj <- get(distance)( t( d )   )
        }else{
          main.dist.obj <- dist( t(d), method=distance )
        }
        attr( main.dist.obj, "method" ) <- distance
      } else stop("unsupported distance specified.")
    } else {
      ## pFeature < 1 or a weightsFeature != NULL
      ## since d is a data matrix, the user wants to sample over the gene features, so main.dist.obj is left as NULL
    }
  }

  for (i in 1:repCount){
    if(verbose){
      message(paste("random subsample",i));
    }
    ## take expression matrix sample, samples and genes
    sample_x = sampleCols( d, pItem, pFeature, weightsItem, weightsFeature )
    this_dist = NA
    if ( ! is.null( main.dist.obj ) ) {
      boot.cols <- sample_x$subcols
      this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
      if ( clusterAlg != "km" ) {
        ## if this isn't kmeans, then convert to a distance object
        this_dist <- as.dist( this_dist )
        attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
      }
    } else {
      ## if main.dist.obj is NULL, then d is a data matrix, and either:
      ##   1) clusterAlg is 'km'
      ##   2) pFeatures < 1 or weightsFeatures have been specified, or
      ##   3) both
      ## so we can't use a main distance object and for every iteration, we will have to re-calculate either
      ##   1) the distance matrix (because we're also sampling the features as well), or
      ##   2) the submat (if using km)

      if ( clusterAlg != "km" )  {
        if ( ! distance %in% acceptable.distance &  ( class(try(get(distance),silent=TRUE))!="function")  ) stop("unsupported distance.")
        if( ( class(try(get(distance),silent=TRUE))=="function") ){
          this_dist <- get(distance)( t( sample_x$submat ) )
        }else{
          if( distance == "pearson" | distance == "spearman"){
            this_dist <- as.dist( 1-cor(sample_x$submat,use=corUse,method=distance) )
          }else{
            this_dist <- dist( t( sample_x$submat ), method= distance  )
          }
        }
        attr( this_dist, "method" ) <- distance
      } else {
        ## if we're not sampling the features, then grab the colslice
        if ( is.null( pFeature ) ||
             ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) {
          this_dist <- d[, sample_x$subcols ]
        } else {
          if ( is.na( sample_x$submat ) ) {
            stop( "error submat is NA" )
          }
          this_dist <- sample_x$submat
        }
      }
    }
    ## cluster samples for HC.
    this_cluster=NA
    if(clusterAlg=="hc"){
      this_cluster = hclust( this_dist, method=innerLinkage)
    }
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix( rep( 1,length(sample_x[[3]])),
                                  mCount,
                                  sample_x[[3]] )
    ##use samples for each k
    for (k in 2:maxK){
      if(verbose){
        message(paste("  k =",k))
      }
      if (i==1){
        ml[[k]] = mConsist #initialize
      }
      this_assignment=NA
      if(clusterAlg=="hc"){
        ##prune to k for hc
        this_assignment = cutree(this_cluster,k)

      }else if(clusterAlg=="kmdist"){
        this_assignment = kmeans(this_dist, k, iter.max = 10, nstart = 1, algorithm = c("Hartigan-Wong") )$cluster

      }else if(clusterAlg=="km"){
        ##this_dist should now be a matrix corresponding to the result from sampleCols
        this_assignment <- kmeans( t( this_dist ),
                                   k,
                                   iter.max = 10,
                                   nstart = 1,
                                   algorithm = c("Hartigan-Wong") )$cluster
      }else if ( clusterAlg == "pam" ) {
        this_assignment <- cluster::pam( x=this_dist,
                                k,
                                diss=TRUE,
                                metric=distance,
                                cluster.only=TRUE )
      } else{
        ##optional cluterArg Hook.
        this_assignment <- get(clusterAlg)(this_dist, k)
      }
      ml[[k]] <- connectivityMatrix( this_assignment,
                                     ml[[k]],
                                     sample_x[[3]] )
    } # END OF INNER K FOR LOOP
  } # END OF OUTER I ITERATION LOOP

  ##consensus fraction
  res = vector(mode="list",maxK)
  for (k in 2:maxK){
    ##fill in other half of matrix for tally and count.
    tmp = triangle(ml[[k]],mode=3)
    tmpCount = triangle(mCount,mode=3)
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
  #message("end fraction")
  return(res)
}

connectivityMatrix <- function( clusterAssignments, m, sampleKey){
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix
  names( clusterAssignments ) <- sampleKey
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId

  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    m <- m + updt
  }
  return(m)
}

sampleCols <- function( d,
                        pSamp=NULL,
                        pRow=NULL,
                        weightsItem=NULL,
                        weightsFeature=NULL ){
  ## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
  ##  if no sampling over the features is performed, the submatrix & sample features are returned as NAs
  ##  to reduce memory overhead
  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  sampleN <- floor(space*pSamp)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem) )
  this_sample <- sampRows <- NA
  if ( inherits( d, "matrix" ) ) {
    if ( (! is.null( pRow ) ) &&
         ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
      ## only sample the rows and generate a sub-matrix if we're sampling over the row/gene/features
      space = nrow(d)
      sampleN = floor(space*pRow)
      sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
      this_sample <- d[sampRows,sampCols]
      dimnames(this_sample) <- NULL
    } else {
      ## do nothing
    }
  }
  return( list( submat=this_sample,
                subrows=sampRows,
                subcols=sampCols ) )
}

COP=function(ml){ # calculates cophenetic distance
  maxK = length(ml)
  copres <- matrix(nrow=(maxK - 1),ncol=2)
  for (ccm in seq(2,maxK,1)){
    x <- ml[[ccm]] # this should be the CC matrix
    d1 <- dist(x)
    hc <- hclust(d1, "average")
    d2 <- cophenetic(hc)
    ans <- cor(d1, d2) # for every cc matrix
    copres[ccm-1,1] <- ccm
    copres[ccm-1,2] <- ans
  }
  copres <- data.frame(copres)
  colnames(copres) <- c('K', 'Cophenetic Coeff')
  return(copres)
}

CDF=function(ml,breaks=100,printres=printres,x1=x1,x2=x2){ # calculates CDF and PAC
  maxK = length(ml) # match with max K
  cdf_res <- matrix(nrow = 10000, ncol = 3)
  i = 1
  for (ccm in seq(2,maxK,1)){
    x <- ml[[ccm]] # this should be the CC matrix
    #x <- unlist(x)
    p <- stats::ecdf(x)
    for (index in seq(0,1,0.01)){
      answer <- p(index)
      cdf_res[i,1] <- index
      cdf_res[i,2] <- answer
      cdf_res[i,3] <- ccm
      i = i + 1
    }
  }
  # CDF plot
  cdf_res2 <- as.data.frame(cdf_res)
  colnames(cdf_res2) <- c('consensusindex', 'CDF', 'k')
  cdf_res2 <- cdf_res2[stats::complete.cases(cdf_res2),]
  p <- ggplot2::ggplot(cdf_res2, ggplot2::aes(x=consensusindex, y=CDF, group=k)) + ggplot2::geom_line(ggplot2::aes(colour = factor(k)), size = 2) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 26, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = 26, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = 26),
          axis.title.y = ggplot2::element_text(size = 26),
          legend.title = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 26),
          plot.title = ggplot2::element_text(size = 26, colour = 'black', hjust = 0.5),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = "Real Data")
  if (printres == TRUE){
    png('CDF.png', height = 14, width = 23, units = 'cm',
         type = 'cairo', res = 900)
  }
  print(p) # print ggplot CDF in main plotting window
  if (printres == TRUE){
    dev.off()
  }
  # Code for calculating and plotting PAC score
  cdf_res3 <- subset(cdf_res2, consensusindex %in% c(x1, x2)) # select the consensus index vals to determine the PAC score
  PAC_res <- matrix(nrow=(maxK - 1),ncol=2)
  j = 2
  r = 1
  for (i in seq(2, nrow(cdf_res3), 2)){
    value1 <- cdf_res3[i,2]
    value2 <- cdf_res3[(i-1),2]
    PAC <- value1-value2
    PAC_res[r,2] <- PAC
    PAC_res[r,1] <- j
    j = j + 1
    r = r + 1
  }
  PAC_res_df <- as.data.frame(PAC_res)
  colnames(PAC_res_df) <- c('K', 'PAC_SCORE')
  # do PAC plot
  p2 <- ggplot2::ggplot(data=PAC_res_df, ggplot2::aes(x=K, y=PAC_SCORE)) + ggplot2::geom_line(colour = "sky blue", size = 2) + ggplot2::geom_point(colour = "black", size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 26, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = 26, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = 26),
          axis.title.y = ggplot2::element_text(size = 26),
          legend.text = ggplot2::element_text(size = 26),
          legend.title = ggplot2::element_text(size = 26),
          plot.title = ggplot2::element_text(size = 26, colour = 'black', hjust = 0.5),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(breaks=c(seq(0,maxK,1))) +
    ggplot2::ylab('PAC Score') +
    ggplot2::xlab('Number of Clusters') +
    ggplot2::labs(title = "Real Data")
  if (printres == TRUE){
    png('PACscore.png', height = 14, width = 20, units = 'cm',
        type = 'cairo', res = 900)
  }
  print(p2)
  if (printres == TRUE){
    dev.off()
  }

  if (printres == TRUE){
    print(p)
    print(p2)
  }

  return(PAC_res_df)
}

myPal = function(n=10){
  #returns n colors
  seq = rev(seq(0,255,by=255/(n)))
  palRGB = cbind(seq,seq,255)
  rgb(palRGB,maxColorValue=255)
}

setClusterColors = function(past_ct,ct,colorU,colorList){
  #description: sets common color of clusters between different K
  newColors = c()
  if(length(colorList)==0){
    #k==2
    newColors = colorU[ct]
    colori=2
  }else{
    newColors = rep(NULL,length(ct))
    colori = colorList[[2]]
    mo=table(past_ct,ct)
    m=mo/apply(mo,1,sum)
    for(tci in 1:ncol(m)){ # for each cluster
      maxC = max(m[,tci])
      pci = which(m[,tci] == maxC)
      if( sum(m[,tci]==maxC)==1 & max(m[pci,])==maxC & sum(m[pci,]==maxC)==1  )  {
        #if new column maximum is unique, same cell is row maximum and is also unique
        ##Note: the greatest of the prior clusters' members are the greatest in a current cluster's members.
        newColors[which(ct==tci)] = unique(colorList[[1]][which(past_ct==pci)]) # one value
      }else{ #add new color.
        colori=colori+1
        newColors[which(ct==tci)] = colorU[colori]
      }
    }
  }
  return(list(newColors,colori,unique(newColors) ))
}

triangle = function(m,mode=1){
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m
  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half
  fm = t(nm)+nm
  diag(fm) = diag(m)
  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]
  if(mode==1){
    return(vm) #vector
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }
}

rmv <- function(n, covmat, rfunc = rnorm, method = c('chol', 'eigen'),...){ # for cholesky method
  m <- nrow(covmat)
  msr <- msr(covmat, method = method)
  x <- matrix(rfunc(n*m, ...), n, m)
  x %*% msr
}

msr <- function(sqrmat, method = c('chol', 'eigen')){ # for cholesky method
  m <- nrow(sqrmat)
  method <- match.arg(method)
  if(method == 'chol'){
    msr <- chol(sqrmat)
  }
  else if(method == 'eigen'){
    ed <- eigen(sqrmat)
    msr <- diag(sqrt(ed$values), m, m) %*% t(ed$vectors)
  }
  else{
    stop('method not recognised')
  }
  return(msr)
}

