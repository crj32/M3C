#' M3C: Monte Carlo Consensus Clustering
#'
#' This function runs M3C, which is a hypothesis testing framework for consensus clustering. The basic
#' idea is to use a multi-core enabled monte carlo simulation to drive the creation of a null distribution
#' of stability scores. The monte carlo simulations maintains the correlation structure of the input data.
#' Then the null distribution is used to compare the reference scores with the real scores
#' and a empirical p value is calculated for every value of K. We also use the relative cluster stability
#' index as an alternative metric which is just based on a comparison against the reference mean, the advantage 
#' being it requires fewer iterations. Small p values are estimated cheaply using a beta distribution that is
#' inferred using parameter estimates from the monte carlo simulation.
#'
#' @param mydata Data frame or matrix: Contains the data, with samples as columns and rows as features
#' @param montecarlo Logical flag: whether to run the monte carlo simulation or not (recommended: TRUE)
#' @param cores Numerical value: how many cores to split the monte carlo simulation over
#' @param iters Numerical value: how many monte carlo iterations to perform (default: 100, recommended: 100-1000)
#' @param maxK Numerical value: the maximum number of clusters to test for, K (default: 10)
#' @param des Data frame: contains annotation data for the input data for automatic reordering (optional)
#' @param ref_method Character string: refers to which reference method to use (recommended: leaving as default)
#' @param repsref Numerical value: how many reps to use for the monte carlo reference data (suggest 100)
#' @param repsreal Numerical value: how many reps to use for the real data (recommended: 100)
#' @param clusteralg String: dictates which algorithm to use for M3C (recommended: leaving as default)
#' @param distance String: dictates which distance metric to use for M3C (recommended: leaving as default)
#' @param pacx1 Numerical value: The 1st x co-ordinate for calculating the pac score from the CDF (default: 0.1)
#' @param pacx2 Numerical value: The 2nd x co-ordinate for calculating the pac score from the CDF (default: 0.9)
#' @param printres Logical flag: whether to print all results into current directory
#' @param printheatmaps Logical flag: whether to print all the heatmaps into current directory
#' @param showheatmaps Logical flag: whether to show the heatmaps on screen (can be slow)
#' @param removeplots Logical flag: whether to remove all plots (recommended: leaving as default)
#' @param seed Numerical value: fixes the seed if you want to repeat results, set the seed to 123 for example here
#'
#' @return A list, containing: 
#' 1) the stability results and 
#' 2) all the output data (another list) 
#' 3) reference stability scores
#' (see vignette for more details on how to easily access)
#' @export
#'
#' @examples
#' res <- M3C(mydata, cores=1, iters=100, ref_method = 'reverse-pca', montecarlo = TRUE,printres = FALSE, 
#' maxK = 10, showheatmaps = FALSE, repsreal = 100, repsref = 100,printheatmaps = FALSE, seed = 123, des = desx)
M3C <- function(mydata, montecarlo = TRUE, cores = 1, iters = 100, maxK = 10,
                              des = NULL, ref_method = c('reverse-pca', 'chol'), repsref = 100, repsreal = 100,
                              clusteralg = c('pam', 'km'), distance = 'euclidean', pacx1 = 0.1, pacx2 = 0.9, printres = FALSE,
                              printheatmaps = FALSE, showheatmaps = FALSE, seed=NULL, removeplots = FALSE){

  if (is.null(seed) == FALSE){
    set.seed(seed)
  }
  
  ref_method <- match.arg(ref_method)
  clusteralg <- match.arg(clusteralg)
  
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
  if (class(mydata) == 'matrix'){
    mydata <- data.frame(mydata)
    colnames(mydata) <- gsub('X', '', colnames(mydata))
  }

  # consensuscluster2 reference or no reference functions

  if (montecarlo == TRUE){
    ## run monte carlo simulation to generate references with same gene-gene correlation structure
    message('running simulations...')
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    invisible(capture.output(pb <- txtProgressBar(min = 0, max = iters, style = 3)))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    ## pre calculations before multi core
    c <- ncol(mydata)
    r <- nrow(mydata)
    if (ref_method == 'reverse-pca'){
      pca1 <- prcomp(t(mydata))
      sds <- colSdColMeans(pca1$x)
    }else if (ref_method == 'chol'){
      if (matrixcalc::is.positive.definite(cov(t(mydata))) == FALSE){
        covm <- Matrix::nearPD(cov(t(mydata)))
        covm <- covm$mat
      }else{
        covm <- cov(t(mydata))
      }
    }
    ## for each loop to use all cores
    ls<-foreach(i = 1:iters, .export=c("ccRun", "CDF", "connectivityMatrix", "M3Cref", "myPal",
                                       "sampleCols", "setClusterColors", "triangle"),
                .packages=c("cluster", "base", "Matrix"), .combine='rbind', .options.snow = opts) %dopar% {

                  if (is.null(seed) == FALSE){
                    set.seed(i)
                  }
                  
                  if (ref_method == 'reverse-pca'){ # reverse PCA method
                    simulated_data <- matrix(rnorm(c*c, mean = 0, sd = sds),nrow = c, ncol = c, byrow = TRUE)
                    null_dataset <- t(t(simulated_data %*% t(pca1$rotation)) + pca1$center) # go back to null distrib
                    mydata2 <- t(as.data.frame(null_dataset))
                  }else if (ref_method == 'chol') { # cholesky method
                    newdata <- matrix(rnorm(c*r), c, r) %*% chol(covm)
                    mydata2<- as.data.frame(t(newdata))
                  }

                  m_matrix <- as.matrix(mydata2)
                  results <- M3Cref(m_matrix,maxK=maxK,reps=repsref,pItem=0.8,pFeature=1,
                                                 clusterAlg=clusteralg, # use pam it is fast
                                                 distance=distance, # with pam always use euclidean
                                                 title = '/home/christopher/Desktop/',
                                                 x1=pacx1, x2=pacx2, printres=FALSE, seed=seed)
                  pacresults <- results$pac_scores$PAC_SCORE
                  return(pacresults)
                }
    close(pb)
    stopCluster(cl)
    message('finished generating reference distribution')
    ## finished monte carlo, results are in ls matrix
    ## real data PAC score calculation
    results2 <- M3Creal(as.matrix(mydata),maxK=maxK,reps=repsreal,pItem=0.8,pFeature=1,
                                     clusterAlg=clusteralg, # use pam it is fast
                                     distance=distance, # with pam always use euclidean
                                     title = '/home/christopher/Desktop/',
                                     printres = printres,
                                     showheatmaps = showheatmaps, printheatmaps = printheatmaps, des = des,
                                     x1=pacx1, x2=pacx2, seed=seed, removeplots=removeplots) # png to file
    real <- results2$pac_scores
    allresults <- results2$allresults

    ## process reference data and calculate scores and p values (simulations are in ls matrix)
    colnames(real)[2] <- 'PAC_REAL'
    real$PAC_REF <- colMeans(ls)
    #real$RCSI <- real$PAC_REF - real$PAC_REAL # alternative rcsi calculation without taking log

    ## if PAC is zero set it to really really small 
    real$PAC_REAL[real$PAC_REAL==0] <- 0.0000001
    pacreal <- real$PAC_REAL
    PACREALLOG <- log(pacreal)
    PACREFLOG <- log(real$PAC_REF)
    real$RCSI <- PACREFLOG - PACREALLOG # calculate RSCI

    ## usual p value derivation
    pvals <- vapply(seq_len(ncol(ls)), function(i) {
      distribution <- as.numeric(ls[,i])
      ((length(distribution[distribution < real$PAC_REAL[i]])) + 1)/(iters+1) # (b+1)/(m+1)=pval
    }, numeric(1))
    real$MONTECARLO_P <- pvals
    
    ## estimate p values using a beta distribution
    variance <- apply(ls, 2, var)
    pvals2 <- vapply(seq_len(nrow(real)), function(i) {
      mean <- real$PAC_REF[i]
      var <- variance[[i]]
      realpac <- real$PAC_REAL[i]
      params2 <- estBetaParams(mu=mean, var=var)
      pbeta(realpac, params2[[1]], params2[[2]])
    }, numeric(1))
    real$BETA_P <- pvals2
    real$P_SCORE <- -log10(real$BETA_P)
    
  if (removeplots == FALSE){ # we are doing the plots
    # plot real vs reference results
    # pac statistic
    px <- ggplot(data=real, aes(x=K, y=RCSI)) + geom_line(colour = "purple", size = 2) + 
      geom_point(colour = "black", size = 3) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 26, colour = 'black'),
            axis.text.x = element_text(size = 26, colour = 'black'),
            axis.title.x = element_text(size = 26),
            axis.title.y = element_text(size = 26),
            legend.text = element_text(size = 26),
            legend.title = element_text(size = 26),
            plot.title = element_text(size = 26, colour = 'black', hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks=c(seq(0,maxK,1))) +
      ylab('RCSI') +
      xlab('K')
    
    # pval score
    col = ifelse(real$P_SCORE > 1.30103,'tomato','black')
    py <- ggplot(data=real, aes(x=K, y=P_SCORE)) + geom_point(colour = col, size = 3) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 26, colour = 'black'),
            axis.text.x = element_text(size = 26, colour = 'black'),
            axis.title.x = element_text(size = 26),
            axis.title.y = element_text(size = 26),
            legend.text = element_text(size = 26),
            legend.title = element_text(size = 26),
            plot.title = element_text(size = 26, colour = 'black', hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks=c(seq(0,maxK,1))) +
      ylab(expression('-log'[10]*'p')) +
      xlab('K') +
      geom_hline(yintercept=1.30103, size=0.75, linetype='dashed', colour='tomato') # 0.05 sig threshold
    
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
  }

  if (montecarlo == FALSE){
    message('running without monte carlo simulations...')
    results2 <- M3Creal(as.matrix(mydata),maxK=maxK,reps=repsreal,pItem=0.8,pFeature=1, 
                                     clusterAlg=clusteralg, # default = pam, others = hc, km
                                     distance=distance, # seed=1262118388.71279,
                                     title = '/home/christopher/Desktop/',
                                     printres = printres, x1=pacx1, x2=pacx2,
                                     showheatmaps = showheatmaps, printheatmaps = printheatmaps, des = des,
                                     seed=seed, removeplots=removeplots) # png to file
    real <- results2$pac_scores
    allresults <- results2$allresults
  }

  if (file.exists('Rplots.pdf') == TRUE){ # random file pops up
    file.remove('Rplots.pdf') # this needs to be removed
    unlink('Rplots.pdf') # this needs to be removed
  }

  if (printres == TRUE){
    write.csv(real, file = 'pacresultfile.csv', row.names = FALSE)
  }

  if (montecarlo == TRUE){
    # return results with monte carlo
    ls <- data.frame(ls)
    row.names(ls) <- gsub('result', 'iteration', row.names(ls))
    colnames(ls) <- c(2:maxK)
    return(list("realdataresults" = allresults, 'scores' = real, 'refpacscores' = ls)) 
  }
  if (montecarlo == FALSE){
    # return results without monte carlo
    return(list("realdataresults" = allresults, 'scores' = real)) 
  }

}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

colSdColMeans <- function(x, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x))
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
  return(sqrt(colVar * n/(n-1)))
}

M3Creal <- function( d=NULL, # function for real data
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
                                  corUse="everything",
                                  showheatmaps=FALSE,
                                  printheatmaps=FALSE,
                                  printres=FALSE,
                                  x1=0.1,
                                  x2=0.9,
                                  des = NULL,
                                  removeplots=removeplots) {
  message('running consensus cluster algorithm for real data...')
  if (is.null(seed) == FALSE){
    set.seed(seed)
  }
  ml <- ccRun(d=d,
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
              corUse=corUse)
  message('finished')
  res=list();
  colorList=list()
  colorM = rbind()
  thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
               "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
               "#bd18ea","#2ef4ca","#f4cced","#f4cc03","#05188a","#e5a25a", 
               "#06f106","#85848f","#000000","#076f25","#93cd7f","#4d0776", 
               "#ffffff")
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
  
  ## loop over each consensus matrix and get the results out
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
    pc=pc[hc$order,]
    pc = rbind(pc,0)
    colcols <- as.factor(as.numeric(as.factor(colorList[[1]])))
    cols <- colorRampPalette(RColorBrewer::brewer.pal(9,'Reds')[1:6])(256)
    if (showheatmaps == TRUE & printheatmaps == FALSE){
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), 
                    col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
    }
    if (showheatmaps == TRUE & printheatmaps == TRUE){
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), 
                    col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      png(paste('K=',tk,'heatmap.png'), height = 12, width = 12, units = 'cm', 
          res = 900, type = 'cairo')
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), 
                    col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      dev.off()
    }
    if (showheatmaps == FALSE & printheatmaps == TRUE){
      png(paste('K=',tk,'heatmap.png'), height = 12, width = 12, units = 'cm', 
          res = 900, type = 'cairo')
      NMF::aheatmap(pc, Colv = as.dendrogram(hc), Rowv = NA, annCol = data.frame(CC = colcols), 
                    col = cols, cexRow = 0, cexCol = 0, annLegend = FALSE)
      dev.off()
    }
    res[[tk]] = list(consensusMatrix=c,consensusTree=hc,consensusClass=ct,ml=ml[[tk]],clrs=colorList)
    colorM = rbind(colorM,colorList[[1]])
    
    ## start code for extracting ordered data out for user (removed the X reformatting)
    usethisorder <- hc$order
    colnames(pc) <- seq(1,ncol(pc))
    pc <- pc[,usethisorder] # get the consensus matrices in the correct order
    m_matrix=as.matrix(d) # this is the input data which we will reorder to match the consensus clusters
    clustering <- as.numeric(as.factor(colorList[[1]]))
    clusteringdf <- data.frame(sample = colnames(m_matrix), cluster = clustering)
    neworder1 <- colnames(m_matrix)[hc$order]
    neworder2 <- gsub('-', '.', neworder1)
    df <- data.frame(m_matrix)
    newdes <- data.frame(ID = colnames(df), consensuscluster = factor(clusteringdf$cluster))
    #colnames(df) <- gsub('X', '', colnames(df))
    #neworder2 <- gsub('X', '', neworder2)
    data <- df[neworder2]
    if (is.null(des) == TRUE){ 
      neworder1 <- gsub('-', '.', neworder1) # cc code is changing - to . so change back
      #vec <- grepl('X', colnames(d)) # check for X's in original colnames if dont exist run this code
      #if (all(vec == FALSE)){ # just changed this to FALSE
        #newdes$ID <- gsub('X', '', newdes$ID) # this creates problem if X in original names
      #}
      newerdes <- newdes[match(neworder1, newdes$ID),]
      annotation <- data.frame(newerdes$consensuscluster)
      row.names(annotation) <- newerdes$ID
      colnames(annotation) <- c('consensuscluster')
    }
    if (is.null(des) == FALSE){
      neworder1 <- gsub('-', '.', neworder1)
      #neworder1 <- gsub('X', '', neworder1)
      des$ID <- gsub('-', '.', des$ID) # this is a problem with formatting of the ids - check out later
      #vec <- grepl('X', colnames(d)) # check for X's in original colnames if dont exist run this code
      #if (all(vec == FALSE)){ # changed to true** then back to false, problems here
      #  newdes$ID <- gsub('X', '', newdes$ID) # this creates problem if X in original names
      #}
      merged <- merge(newdes, des, by = 'ID') # name to merge by
      newdes <- merged
      newerdes <- newdes[match(neworder1, newdes$ID),]
      annotation <- newerdes
      row.names(annotation) <- newerdes$ID
      annotation$ID <- NULL
    }
    
    # change the numbering of the consensus clusters so they make sense
    vec <- annotation$consensuscluster
    maximum <- length(levels(vec))
    seq <- seq(1,maximum)
    j = 1
    vec <- as.character(vec)
    vec2 <- c()
    for (i in seq(1,length(vec))){
      if (i == 1){ # if we are at the first element
        vec2[i] <- 1
      }
      if (i > 1){ # if we are at any other element
        x <- vec[i]
        y <- vec[i-1]
        if (x == y){
          vec2[i] <- seq[j]
        }
        if (x != y){
          j = j + 1
          vec2[i] <- seq[j]
        }
      }
    }
    annotation$consensuscluster <- as.factor(vec2)
    
    newList <- list("consensus_matrix" = pc, 'ordered_data' = data, 'ordered_annotation' = annotation) # you can remove ml
    resultslist[[tk]] <- newList
  }
  pac_res <- CDF(ml, printres=printres, x1=x1, x2=x2, removeplots=removeplots) # this runs the new CDF function with PAC score
  res[[1]] = colorM
  listxyx <- list("allresults" = resultslist, 'pac_scores' = pac_res) # use a name list, one item is a list of results
  return(listxyx)
  print('finished this function')
}

M3Cref <- function( d=NULL, # function for reference data
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
  message('running consensus cluster algorithm for reference data...') # this is the main function that takes the vast majority of the time
  ml <- ccRun(d=d,
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
              corUse=corUse)
  message('finished.')
  pac_res <- CDF(ml, printres=FALSE, x1=x1, x2=x2, removeplots=TRUE) # this runs the new CDF function with PAC score
  newList <- list('pac_scores' = pac_res) # now returning a fairly coherant list of results
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
                   corUse=NULL) {
  
  ## setting up initial objects
  m = vector(mode='list', repCount)
  ml = vector(mode="list",maxK)
  n <- ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  mCount = mConsist = matrix(c(0),ncol=n,nrow=n)
  ml[[1]] = c(0);
  acceptable.distance <- c( "euclidean")
  main.dist.obj <- NULL
  
  ## if pam you need to make the distance matrix first
  if ( clusterAlg == "pam" ){
    main.dist.obj <- dist( t(d), method=distance )
  }
  
  ## start the resampling loop
  for (i in 1:repCount){ 
    ## sample the input data matrix
    sample_x = sampleCols(d, pItem, pFeature, weightsItem, weightsFeature)
    this_dist = NA
    if (clusterAlg == "pam"){ # if algorithm equals PAM do this
      boot.cols <- sample_x$subcols
      this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
      this_dist <- as.dist( this_dist )
      attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
    }else if (clusterAlg == 'km') { # algorithm equals KMEANS then do this
      this_dist <- d[, sample_x$subcols ]
    }
    ## mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix( rep( 1,length(sample_x[[3]])),
                                  mCount,
                                  sample_x[[3]] )
    ## loop over different values of K
    for (k in 2:maxK){
      if (i==1){
        ml[[k]] = mConsist
      }
      this_assignment=NA
      if(clusterAlg=="km"){
        this_assignment <- kmeans(t(this_dist),k,iter.max = 10,nstart = 1,
                                  algorithm = c("Hartigan-Wong") )$cluster
      }else if ( clusterAlg == "pam" ) {
        this_assignment <- cluster::pam(x=this_dist,k,diss=TRUE,metric=distance,cluster.only=TRUE)
      } 
      ml[[k]] <- connectivityMatrix(this_assignment,ml[[k]],sample_x[[3]]) 
    }
  }
  
  ##consensus fraction
  res = vector(mode="list",maxK)
  for (k in 2:maxK){
    tmp = triangle(ml[[k]],mode=3)
    tmpCount = triangle(mCount,mode=3)
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
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
  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  sampleN <- floor(space*pSamp)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem) )
  this_sample <- sampRows <- NA
  if ( inherits( d, "matrix" ) ) {
    if ( (! is.null( pRow ) ) &&
         ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
      space = nrow(d)
      sampleN = floor(space*pRow)
      sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
      this_sample <- d[sampRows,sampCols]
      dimnames(this_sample) <- NULL
    } else {
      #
    }
  }
  return( list( submat=this_sample,
                subrows=sampRows,
                subcols=sampCols ) )
}

CDF=function(ml,breaks=100,printres=printres,x1=x1,x2=x2,
             removeplots=removeplots){ # calculate CDF and PAC score
  
  ## calculate CDF
  maxK = length(ml) # match with max K
  cdf_res <- matrix(nrow = 10000, ncol = 3)
  i = 1
  for (ccm in seq(2,maxK,1)){ 
    x <- ml[[ccm]] # this should be the CC matrix
    x <- x[lower.tri(x)]
    p <- ecdf(x)
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
  cdf_res2 <- cdf_res2[complete.cases(cdf_res2),]
  
  if (removeplots==FALSE){
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
  }

  ## vectorised PAC score calculation
  cdf_res3 <- subset(cdf_res2, consensusindex %in% c(x1, x2)) # select the consensus index vals to determine the PAC score
  value1 <- cdf_res3[seq(2, nrow(cdf_res3), 2), 2]
  value2 <- cdf_res3[seq(1, nrow(cdf_res3), 2), 2]
  PAC <- value1 - value2
  PAC_res_df <- data.frame(K=seq(2, maxK), PAC_SCORE=PAC)
  
  if (removeplots==FALSE){
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
  }
  
  return(PAC_res_df)
}

myPal = function(n=10){
  seq = rev(seq(0,255,by=255/(n)))
  palRGB = cbind(seq,seq,255)
  rgb(palRGB,maxColorValue=255)
}

setClusterColors = function(past_ct,ct,colorU,colorList){
  newColors = c()
  if(length(colorList)==0){
    newColors = colorU[ct]
    colori=2
  }else{
    newColors = rep(NULL,length(ct))
    colori = colorList[[2]]
    mo=table(past_ct,ct)
    m=mo/apply(mo,1,sum)
    for(tci in 1:ncol(m)){
      maxC = max(m[,tci])
      pci = which(m[,tci] == maxC)				
      if( sum(m[,tci]==maxC)==1 & max(m[pci,])==maxC & sum(m[pci,]==maxC)==1  )  {
        newColors[which(ct==tci)] = unique(colorList[[1]][which(past_ct==pci)])
      }else{
        colori=colori+1
        newColors[which(ct==tci)] = colorU[colori]
      }
    }
  }
  return(list(newColors,colori,unique(newColors) ))
}

triangle = function(m,mode=1){
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
