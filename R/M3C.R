#' M3C: Monte Carlo Reference-based Consensus Clustering
#'
#' This is the M3C core function, which is a reference-based consensus clustering algorithm. The basic
#' idea is to use a multi-core enabled Monte Carlo simulation to drive the creation of a null distribution
#' of stability scores. The Monte Carlo simulations maintains the feature correlation structure of the 
#' input data. Then the null distribution is used to compare the reference scores with the real scores
#' and an empirical p value is calculated for every value of K to test the null hypothesis K=1. We derive 
#' the Relative Cluster Stability Index (RCSI) as a metric for selecting K, which is based on a 
#' comparison against the reference mean. A fast alternative is also included that includes a penalty
#' term to prevent overestimation of K, we call regularised consensus clustering.
#'
#' @param mydata Data frame or matrix: Contains the data, with samples as columns and rows as features
#' @param cores Numerical value: how many cores to split the monte carlo simulation over
#' @param iters Numerical value: how many Monte Carlo iterations to perform (default: 25, recommended: 5-100)
#' @param maxK Numerical value: the maximum number of clusters to test for, K (default: 10)
#' @param des Data frame: contains annotation data for the input data for automatic reordering
#' @param ref_method Character string: refers to which reference method to use
#' @param repsref Numerical value: how many resampling reps to use for reference (default: 100, recommended: 100-250)
#' @param repsreal Numerical value: how many resampling reps to use for real data (default: 100, recommended: 100-250)
#' @param clusteralg String: dictates which inner clustering algorithm to use (default: PAM)
#' @param pacx1 Numerical value: The 1st x co-ordinate for calculating the pac score from the CDF (default: 0.1)
#' @param pacx2 Numerical value: The 2nd x co-ordinate for calculating the pac score from the CDF (default: 0.9)
#' @param removeplots Logical flag: whether to remove all plots from view
#' @param fsize Numerical value: determines the font size of the ggplot2 plots
#' @param lthick Numerical value: determines the line thickness of the ggplot2 plot
#' @param dotsize Numerical value: determines the dotsize of the ggplot2 plot
#' @param method Numerical value: 1 refers to the Monte Carlo simulation method, 2 to regularised consensus clustering
#' @param seed Numerical value: specifies seed, set to NULL for different results each time
#' @param tunelambda Logical flag: whether to tune lambda or not
#' @param lseq Numerical vector: vector of lambda values to tune over (default = seq(0.05,0.1,by=0.01))
#' @param lambdadefault Numerical value: if not tuning fixes the default (default: 0.1)
#' @param silent Logical flag: whether to remove messages or not
#' @param objective Character string: whether to use 'PAC' or 'entropy' objective function (default = entropy)
#' @param pItem Numerical value: the fraction of points to resample each iteration (default: 0.8)
#'
#' @return A list, containing: 
#' 1) the stability results and 
#' 2) all the output data (another list) 
#' 3) reference stability scores
#' (see vignette for more details on how to easily access)
#' @export
#'
#' @examples
#' res <- M3C(mydata)
M3C <- function(mydata, cores = 1, iters = 25, maxK = 10, pItem = 0.8,
                des = NULL, ref_method = c('reverse-pca', 'chol'), repsref = 100, repsreal = 100,
                clusteralg = c('pam', 'km', 'spectral', 'hc'), pacx1 = 0.1, 
                pacx2 = 0.9, seed=123, objective='entropy', removeplots = FALSE,
                silent = FALSE, fsize = 18, method = 1, lambdadefault = 0.1, tunelambda = TRUE,
                lseq = seq(0.02,0.1,by=0.02), lthick=2, dotsize=3){
  
  if (is.null(seed) == FALSE){
    set.seed(seed)
  }
  
  ref_method <- match.arg(ref_method)
  clusteralg <- match.arg(clusteralg)
  distance <- 'euclidean' # always use this
  
  if (method == 1){
    montecarlo <- TRUE 
  }else if (method == 2){
    montecarlo <- FALSE 
  }
  
  if (silent != TRUE){
    message('***M3C***')
    if (method == 1){
      message('method: Monte Carlo simulation')
    }else if (method == 2){
      message('method: regularised consensus clustering')
    }
    if (objective == 'entropy'){
      message('objective: entropy')
    }else if (objective == 'PAC'){
      message('objective: pac')
    }
    message(paste('clustering algorithm:',clusteralg))
  }
  
  # error handling of input variables
  
  if ( inherits( mydata,"ExpressionSet" ) ) {
    mydata <- exprs(mydata)
  }
  if (method == 1){
    if (ncol(mydata) > nrow(mydata)){
      if (silent != TRUE){
        message('samples(columns) exceeds features(rows), switching to Cholesky decomposition reference')
      }
      ref_method = 'chol' # this works when variables < samples
    }
  }
  if (is.null(des) == TRUE){ # user does not give extra annotation data
    if (silent != TRUE){
      message('annotation: none')
    }
  }
  if (is.null(des) == FALSE){ # user does give extra annotation data
    if (silent != TRUE){
      message('annotation: yes')
    }
    '%!in%' <- function(x,y)!('%in%'(x,y))
    if("ID" %!in% colnames(des))
    {
      stop('in the supplied annotation object for reordering, ther is no \'ID\' column')
    }
    if (all(colnames(mydata) %in% des$ID) == FALSE){
      stop('the IDs in the annotation do not match column IDs in data')
    }
    if (nrow(des) != ncol(mydata)){
      stop('the dimensions of your annotation object do not match data object')
    }
  }
  if (is(mydata,"matrix")){
    mydata <- data.frame(mydata)
    colnames(mydata) <- gsub('X', '', colnames(mydata))
  }
  
  # M3C reference or no reference functions
  
  if (montecarlo == TRUE){
    ## run monte carlo simulation to generate references with same gene-gene correlation structure
    if (silent != TRUE){
      message('running simulations...')
    }
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
        done <- cov.shrink(t(mydata), verbose=FALSE)
        attr(done, "lambda") <- NULL
        attr(done, "lambda.estimated") <- NULL
        attr(done, "class") <- NULL
        attr(done, "lambda.var") <- NULL
        attr(done, "lambda.var.estimated") <- NULL
        covm <- as.matrix(done)
      }else{
        covm <- cov(t(mydata))
      }
    }
    ## for each loop to use all cores
    ls<-foreach(i = 1:iters, .export=c("ccRun", "CDF", "connectivityMatrix", "M3Cref","entropy",
                                       "sampleCols", "triangle", "rbfkernel"),
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
                  results <- M3Cref(m_matrix,maxK=maxK,reps=repsref,pItem=pItem,pFeature=1,
                                    clusterAlg=clusteralg, # use pam it is fast
                                    distance=distance, # with pam always use euclidean
                                    title = '/home/christopher/Desktop/',
                                    x1=pacx1, x2=pacx2, seed=seed,
                                    silent=silent, objective=objective)
                  pacresults <- results$pac_scores$PAC_SCORE
                  return(pacresults)
                }
    close(pb)
    stopCluster(cl)
    if (silent != TRUE){
      message('done.')
    }
    ## finished monte carlo, results are in ls matrix
    ## real data PAC score calculation
    results2 <- M3Creal(as.matrix(mydata),maxK=maxK,reps=repsreal,pItem=pItem,pFeature=1,
                        clusterAlg=clusteralg, # use pam it is fast
                        distance=distance, # with pam always use euclidean
                        title = '/home/christopher/Desktop/',
                        des = des, lthick=lthick, dotsize=dotsize,
                        x1=pacx1, x2=pacx2, seed=seed, removeplots=removeplots, silent=silent,
                        fsize=fsize,method=method, objective=objective) # png to file
    real <- results2$pac_scores
    allresults <- results2$allresults
    plots <- results2$plots
    
    ## process reference data and calculate scores and p values (simulations are in ls matrix)
    colnames(real)[2] <- 'PAC_REAL'
    real$PAC_REF <- colMeans(ls)
    
    ## if PAC/entropy is zero set it to really really small 
    ptemp <- real$PAC_REAL
    ptemp[ptemp==0] <- 0.0001 ## changed
    pacreal <- ptemp
    
    ## old RCSI calculation log of mean 
    # PACREALLOG <- log(pacreal)
    # PACREFLOG <- log(real$PAC_REF)
    # real$RCSI <- PACREFLOG - PACREALLOG # calculate RSCI
    
    diffM <- sweep(log(ls),2,log(pacreal))
    #diffM <- sweep(ls,2,pacreal)
    real$RCSI <- colMeans(diffM)
    real$RCSI_SE <- (apply(diffM, 2, sd))/sqrt(nrow(ls))
    
    ## usual p value derivation
    pvals <- vapply(seq_len(ncol(ls)), function(i) {
      distribution <- as.numeric(ls[,i])
      ((length(distribution[distribution < real$PAC_REAL[i]])) + 1)/(iters+1) # (b+1)/(m+1)=pval
    }, numeric(1))
    real$MONTECARLO_P <- pvals
    
    if (objective == 'PAC'){
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
    }else if (objective == 'entropy'){
      ## estimate p values using a beta distribution
      variance <- apply(ls, 2, sd)
      pvals2 <- vapply(seq_len(nrow(real)), function(i) {
        mean <- real$PAC_REF[i]
        var <- variance[[i]]
        realpac <- real$PAC_REAL[i]
        pnorm(realpac, mean=mean,sd=var)
      }, numeric(1))
      real$NORM_P <- pvals2
      real$P_SCORE <- -log10(real$NORM_P)
      # fix names
      colnames(real)[2:3]<-c('ENTROPY_REAL','ENTROPY_REF')
    }

    # plot real vs reference results
    # RCSI
    px <- ggplot(data=real, aes(x=K, y=RCSI)) + geom_line(colour = "slateblue", size = lthick) + 
      geom_point(colour = "slateblue", size = dotsize) +
      theme_bw() +
      theme(axis.text.y = element_text(size = fsize, colour = 'black'),
            axis.text.x = element_text(size = fsize, colour = 'black'),
            axis.title.x = element_text(size = fsize),
            axis.title.y = element_text(size = fsize),
            legend.text = element_text(size = fsize),
            legend.title = element_text(size = fsize),
            plot.title = element_text(size = fsize, colour = 'black', hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks=c(seq(0,maxK,1))) +
      ylab('RCSI') +
      xlab('K') + 
      geom_errorbar(aes(ymin=RCSI-qnorm(0.975)*RCSI_SE, ymax=RCSI+qnorm(0.975)*RCSI_SE), 
                    width=.1, colour = "slateblue",size=lthick)
    # pval score
    col = ifelse(real$P_SCORE > 1.30103,'tomato','black')
    py <- ggplot(data=real, aes(x=K, y=P_SCORE)) + geom_point(colour = col, size = dotsize) +
      theme_bw() +
      theme(axis.text.y = element_text(size = fsize, colour = 'black'),
            axis.text.x = element_text(size = fsize, colour = 'black'),
            axis.title.x = element_text(size = fsize),
            axis.title.y = element_text(size = fsize),
            legend.text = element_text(size = fsize),
            legend.title = element_text(size = fsize),
            plot.title = element_text(size = fsize, colour = 'black', hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks=c(seq(0,maxK,1))) +
      ylab(expression('-log'[10]*'p')) +
      xlab('K') #+
    #geom_hline(yintercept=1.30103, size=0.75, linetype='dashed', colour='tomato') # 0.05 sig threshold
    if (removeplots){
      # not much to do here
    }else{
      print(py)
      print(px)
    }
    plots[[3]] <- py
    plots[[4]] <- px
  }
  
  if (montecarlo == FALSE){
    results2 <- M3Creal(as.matrix(mydata),maxK=maxK,reps=repsreal,pItem=pItem,pFeature=1, 
                        clusterAlg=clusteralg, # default = pam, others = hc, km
                        distance=distance, # seed=1262118388.71279,
                        title = '/home/christopher/Desktop/',
                        x1=pacx1, x2=pacx2,
                        des = des,lthick=lthick, dotsize=dotsize,
                        seed=seed, removeplots=removeplots, 
                        silent=silent, 
                        method=method, fsize=fsize, lambda=lambda, objective=objective) # png to file
    real <- results2$pac_scores
    allresults <- results2$allresults
    plots <- results2$plots
  }
  
  if (file.exists('Rplots.pdf') == TRUE){ # random file pops up
    file.remove('Rplots.pdf') # this needs to be removed
    unlink('Rplots.pdf') # this needs to be removed
  }
  
  ### tuning lambda
  if (method == 2 & tunelambda == TRUE){
    llv <- c()
    i = 1
    #lseq <- seq(0.05,0.1,by=0.01)
    ## diff matrix
    for (lambda in lseq){
      if (silent == FALSE){
        message(paste('tuning lambda:',lambda))
      }
      ent <- log(real$PAC_SCORE)+lambda*real$K
      min <- which.min(ent)+1
      #message(paste('K:',min))
      ll <- getl(allresults,min,clusteralg=clusteralg)
      if (silent == FALSE){
        message(paste('K =',min,', log-likelihood',ll))
      }
      llv[i] <- ll
      i = i + 1
    }
    if (silent == FALSE){
      message(paste('optimal lambda:',lseq[tail(which(llv==max(llv)),1)])) # last max
    }
    lambdadefault <- lseq[tail(which(llv==max(llv)),1)]
  }
  
  ### penalised method
  if (objective == 'entropy' & method == 2){
    colnames(real)[2] <- 'ENTROPY'
    real$PCSI <- log(real$ENTROPY)+lambdadefault*real$K
    pz <- PCSI_plot(real,fsize=fsize,maxK=maxK,lthick=lthick, dotsize=dotsize)
    if (removeplots == FALSE){
      print(pz)
    }
    plots[[3]] <- pz
  }else if (objective == 'PAC' & method == 2){
    real$PCSI <- log(real$PAC_SCORE)+lambdadefault*real$K
    pz2 <- PCSI_plot(real,fsize=fsize,maxK=maxK,lthick=lthick, dotsize=dotsize)
    if (removeplots == FALSE){
      print(pz2)
    }
    plots[[3]] <- pz2
  }
  
  ### ground truth compare method
  if (method == 1){
    optk <- which.max(real$RCSI)+1
  }else if (method == 2){
    optk <- which.min(real$PCSI)+1
  }
  
  if (silent != TRUE){
    # print optimal K
    if (method == 1){
      message(paste('optimal K:',optk))
    }else if (method == 2){
      message(paste('optimal K:',optk))
    }
  }
  # get optk assignments out
  assignments <- as.numeric(allresults[[optk]]$assignments)
  
  if (montecarlo == TRUE){
    # return results with monte carlo
    ls <- data.frame(ls)
    row.names(ls) <- gsub('result', 'iteration', row.names(ls))
    colnames(ls) <- c(2:maxK)
    return(list("realdataresults" = allresults, 'scores' = real, 'refpacscores' = ls, 
                'assignments' = assignments, 'plots'=plots))
  }
  if (montecarlo == FALSE){
    # return results without monte carlo
    return(list("realdataresults" = allresults, 'scores' = real, 'assignments'= assignments,
                "plots"=plots)) 
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
                     x1=0.1,
                     x2=0.9,lthick=lthick, dotsize=dotsize,
                     des = NULL,
                     removeplots=removeplots,
                     silent=silent,
                     fsize=fsize,objective=objective,
                     method=method,
                     lambda=lambda) {
  if (silent != TRUE){
    message('running consensus cluster algorithm for real data...')
  }
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
  if (silent != TRUE){
    message('done.')
  }
  ## loop over each consensus matrix and get the results out
  resultslist <- list() # holds all results for each value of K
  for (tk in 2:maxK){
    fm = ml[[tk]]
    hc = hclust( as.dist( 1 - fm ), method=finalLinkage)
    ct = cutree(hc,tk)
    names(ct) = colnames(d)
    if(is(d,"dist")){
      names(ct) = colnames(as.matrix(d))
    }
    c = fm
    pc = c
    pc = pc[hc$order,]
    colcols <- as.factor(as.numeric(ct))
    #pc = rbind(pc,0)
    ## start code for extracting ordered data out for user (removed the X reformatting)
    usethisorder <- hc$order
    colnames(pc) <- seq(1,ncol(pc))
    pc <- pc[,usethisorder] # get the consensus matrices in the correct order
    m_matrix=as.matrix(d) # this is the input data which we will reorder to match the consensus clusters
    clustering <- colcols
    clusteringdf <- data.frame(sample = colnames(m_matrix), cluster = clustering)
    neworder1 <- colnames(m_matrix)[hc$order]
    neworder2 <- gsub('-', '.', neworder1)
    df <- data.frame(m_matrix)
    newdes <- data.frame(ID = colnames(df), consensuscluster = factor(clusteringdf$cluster))
    data <- df[neworder2]
    if (is.null(des) == TRUE){ 
      neworder1 <- gsub('-', '.', neworder1) 
      newerdes <- newdes[match(neworder1, newdes$ID),]
      annotation <- data.frame(newerdes$consensuscluster)
      row.names(annotation) <- newerdes$ID
      colnames(annotation) <- c('consensuscluster')
    }
    if (is.null(des) == FALSE){
      neworder1 <- gsub('-', '.', neworder1)
      des$ID <- gsub('-', '.', des$ID) # this is a problem with formatting of the ids - check out later
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
    
    newList <- list("consensus_matrix" = pc, 'ordered_data' = data, 'ordered_annotation' = annotation, 
                    "assignments" = ct) # you can remove ml
    resultslist[[tk]] <- newList
  }
  
  pac_res <- CDF(ml, x1=x1, x2=x2, removeplots=removeplots, fsize=fsize, lthick=lthick, dotsize=dotsize,
                 method=method, lambda=lambda,objective=objective,returnplots=TRUE) # this runs the new CDF function with PAC score
  
  listxyx <- list("allresults" = resultslist, 'pac_scores' = pac_res$data, 
                  'plots' = pac_res$plots) # use a name list, one item is a list of results
  return(listxyx)
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
                    seed=NULL,
                    silent=silent,objective=objective) {
  if (is.null(seed) == FALSE){
    set.seed(seed)
  }
  if (silent != TRUE){
    message('running consensus cluster algorithm for reference data...') # this is the main function that takes the vast majority of the time
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
  if (silent != TRUE){
    message('done.')
  }
  pac_res <- CDF(ml, x1=x1, x2=x2, removeplots=TRUE, fsize=18, method=1,lthick=1, dotsize=1,
                 objective=objective, returnplots=FALSE) # this runs the new CDF function with PAC score
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
  
  if ( clusterAlg == "pam" ){ # if pam you need to make the distance matrix first
    main.dist.obj <- dist( t(d), method=distance )
  }else if (clusterAlg == 'spectral'){ # do affinity matrix calculation and resample that
    affinitymatrixraw <- rbfkernel(d)
    colnames(affinitymatrixraw) <- colnames(d) # note the order may have changed here
    rownames(affinitymatrixraw) <- colnames(d)
  }else if (clusterAlg == 'hc'){
    main.dist.obj <- dist(t(d),method=distance )
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
    }else if (clusterAlg == 'km') { # if algorithm equals KMEANS then do this
      this_dist <- d[, sample_x$subcols ]
    }else if (clusterAlg == 'spectral'){ # if algorithm equals SPECTRAL do this
      affinitymatrix <- affinitymatrixraw[sample_x$subcols , sample_x$subcols ] # sample beforehand
    }else if (clusterAlg == 'hc'){
      boot.cols <- sample_x$subcols
      this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
      this_dist <- as.dist( this_dist )
      attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
    }
    
    ## mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix( rep( 1,length(sample_x[[3]])),
                                  mCount,
                                  sample_x[[3]] )
    
    ## loop over different values of K
    for (k in 2:maxK){
      #print(k)
      if (i==1){
        ml[[k]] = mConsist
      }
      this_assignment=NA
      if(clusterAlg=="km"){
        this_assignment <- kmeans(t(this_dist),k,iter.max = 10,nstart = 1,
                                  algorithm = c("Hartigan-Wong") )$cluster
      }else if ( clusterAlg == "pam" ) {
        this_assignment <- cluster::pam(x=this_dist,k,diss=TRUE,metric=distance,cluster.only=TRUE)
        #print(this_assignment)
      }else if ( clusterAlg == 'spectral'){
        centers <- k
        m <- nrow(affinitymatrix)
        nc <-  centers
        dv <- 1/sqrt(rowSums(affinitymatrix))
        l <- dv * affinitymatrix %*% diag(dv)
        xi <- eigen(l)$vectors[,1:nc]
        yi <- xi/sqrt(rowSums(xi^2))
        if (any(is.na(yi))){ # fill in columns with column mean when NA
          for(iii in 1:ncol(yi)){
            yi[is.na(yi[,iii]), iii] <- mean(yi[,iii], na.rm = TRUE)
          }
        }
        res <- NULL
        while( is.null(res) ) { # this will hang if there is a problem
          try(
            res <- kmeans(yi, centers)
          )
        } 
        #
        this_assignment <- res$cluster
        names(this_assignment) <- colnames(affinitymatrix)
      }else if (clusterAlg == 'hc'){
        this_cluster <- hclust( this_dist, method='ward.D2')
        this_assignment <- cutree(this_cluster,k)
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

CDF=function(ml,breaks=100,x1=x1,x2=x2,lthick=lthick, dotsize=dotsize,
             removeplots=removeplots,fsize=18,method=1,
             lambda=0.1,objective=objective,returnplots=FALSE){ # calculate CDF and PAC score
  
  plist <- list() # list of plots to save
  
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
  
  if (returnplots == TRUE){
    p <- ggplot2::ggplot(cdf_res2, ggplot2::aes(x=consensusindex, y=CDF, group=k)) + ggplot2::geom_line(ggplot2::aes(colour = factor(k)), size = lthick) + ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = fsize, colour = 'black'),
                     axis.text.x = ggplot2::element_text(size = fsize, colour = 'black'),
                     axis.title.x = ggplot2::element_text(size = fsize),
                     axis.title.y = ggplot2::element_text(size = fsize),
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = fsize),
                     plot.title = ggplot2::element_text(size = fsize, colour = 'black', hjust = 0.5),
                     panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(title = "Real Data")
    if (removeplots == FALSE){
      print(p)
    }
    plist[[1]] <- p
  }
  
  if (objective == 'PAC'){
    ## vectorised PAC score calculation
    cdf_res3 <- subset(cdf_res2, consensusindex %in% c(x1, x2)) # select the consensus index vals to determine the PAC score
    value1 <- cdf_res3[seq(2, nrow(cdf_res3), 2), 2]
    value2 <- cdf_res3[seq(1, nrow(cdf_res3), 2), 2]
    PAC <- value1 - value2
  }
  
  if (objective == 'PAC'){
    ccol <- 'skyblue'
    llab <- 'PAC'
    ## code for penalised PAC
    if (method == 2){
      ## make results data frame
      PAC_res_df <- data.frame(K=seq(2, maxK), PAC_SCORE=PAC)
    }else if (method == 1){ # same as before
      ## make results data frame
      PAC_res_df <- data.frame(K=seq(2, maxK), PAC_SCORE=PAC)
    }
  }else{
    # entropy function
    S <- c()
    for (ccm in seq(2,maxK)){
      cmm <- ml[[ccm]]
      S <- c(S,entropy(cmm))
    }
    S <- data.frame(K=seq(2, maxK),PAC_SCORE=S)
    PAC_res_df <- S
    ccol <- 'black'
    llab <- 'Entropy'
  }
  
  if (returnplots == TRUE){
    p2 <- ggplot2::ggplot(data=PAC_res_df, ggplot2::aes(x=K, y=PAC_SCORE)) + ggplot2::geom_line(colour = ccol, size = lthick) + ggplot2::geom_point(colour = "black", size = dotsize) +
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = fsize, colour = 'black'),
                     axis.text.x = ggplot2::element_text(size = fsize, colour = 'black'),
                     axis.title.x = ggplot2::element_text(size = fsize),
                     axis.title.y = ggplot2::element_text(size = fsize),
                     legend.text = ggplot2::element_text(size = fsize),
                     legend.title = ggplot2::element_text(size = fsize),
                     plot.title = ggplot2::element_text(size = fsize, colour = 'black', hjust = 0.5),
                     panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_x_continuous(breaks=c(seq(0,maxK,1))) +
      ggplot2::ylab(llab) +
      ggplot2::xlab('K') +
      ggplot2::labs(title = "Real Data")
    if (removeplots == FALSE){
      print(p2)
    }
    plist[[2]] <- p2
  }
  
  if (returnplots == TRUE){
    cdfl <- list(data=PAC_res_df,plots=plist)
  }else if (returnplots == FALSE){
    cdfl <- PAC_res_df
  }
  
  return(cdfl)
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

clust.medoid = function(i, distmat, clusters) {
  ind = (clusters == i)
  names(which.min(rowSums( distmat[ind, ind, drop = FALSE] )))
}

rbfkernel <- function (mat, K=3) { # calculate gaussian kernel with local sigma
  n <- ncol(mat)
  NN <- K # nearest neighbours (2-3)
  dm <- as.matrix(dist(t(mat)))
  kn <- c() # find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    sortedvec <- as.numeric(sort.int(dm[i, ]))
    sortedvec <- sortedvec[!sortedvec == 0]
    kn <- c(kn, sortedvec[NN])
  }
  sigmamatrix <- kn %o% kn # make the symmetrical matrix of kth nearest neighbours distances
  upper <- -dm^2 # calculate the numerator beforehand
  out <- matrix(nrow = n, ncol = n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      lowerval <- sigmamatrix[i, j] # retrieve sigma
      upperval <- upper[i, j]
      out[i, j] <- exp(upperval / (lowerval)) # calculate local affinity
    }
  }
  # reflect, diag = 1
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  ## return kernel
  return(out)
}

entropy <- function(m){
  # Zhao, Feng, et al. "Spectral clustering with eigenvector selection 
  # based on entropy ranking." Neurocomputing 73.10-12 (2010): 1704-1717.
  m[m==1] <- NA
  m[m==0] <- NA
  m <- m[upper.tri(m)]
  N <- length(m[!is.na(m)])
  s <- -sum(m*log(m)+(1-m)*log(1-m),na.rm=TRUE)
  if (is.na(s)){
    s<-0
  }
  return(s)
}

PCSI_plot <- function(real,fsize=fsize,maxK=maxK,lthick=lthick, dotsize=dotsize){
  p3 <- ggplot2::ggplot(data=real, ggplot2::aes(x=K, y=PCSI)) + ggplot2::geom_line(colour = "slateblue", size = lthick) + ggplot2::geom_point(colour = "black", size = dotsize) +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fsize, colour = 'black'),
                   axis.text.x = ggplot2::element_text(size = fsize, colour = 'black'),
                   axis.title.x = ggplot2::element_text(size = fsize),
                   axis.title.y = ggplot2::element_text(size = fsize),
                   legend.text = ggplot2::element_text(size = fsize),
                   legend.title = ggplot2::element_text(size = fsize),
                   plot.title = ggplot2::element_text(size = fsize, colour = 'black', hjust = 0.5),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(breaks=c(seq(0,maxK,1))) +
    ggplot2::ylab('PCSI') +
    ggplot2::xlab('K') +
    ggplot2::labs(title = "Real Data")
}

likelihood <- function(m,m2){
  m <- m[upper.tri(m)]
  m2 <- m2[upper.tri(m2)]
  ma <- m[m==1]
  m2a <- m2[m==1]
  L <- sum(log(m2a*ma+(1-ma)*(1-m2a))) # m = ground truth probs, m2 = perturbed probs
  L <- L/length(ma)
  return(L)
}

getl <- function(allresults,k,clusteralg=clusteralg){
  ## get ground truth
  if (clusteralg == 'pam'){
    clust <- cluster::pam(t(allresults[[k]]$ordered_data),k=k)
    clus <- clust$clustering
  }else if (clusteralg == 'km'){
    clus <- kmeans(t(allresults[[k]]$ordered_data),k,iter.max = 10,nstart = 1,
                              algorithm = c("Hartigan-Wong") )$cluster
  }else if (clusteralg == 'hc'){
    this_dist <- dist(t(allresults[[k]]$ordered_data))
    this_cluster <- hclust( this_dist, method='ward.D2')
    clus <- cutree(this_cluster,k)
  }else if (clusteralg == 'spectral'){
    affinitymatrix <- rbfkernel(allresults[[k]]$ordered_data)
    colnames(affinitymatrix) <- colnames(allresults[[k]]$ordered_data) # note the order may have changed here
    rownames(affinitymatrix) <- colnames(allresults[[k]]$ordered_data)
    centers <- k
    m <- nrow(affinitymatrix)
    nc <-  centers
    dv <- 1/sqrt(rowSums(affinitymatrix))
    l <- dv * affinitymatrix %*% diag(dv)
    xi <- eigen(l)$vectors[,1:nc]
    yi <- xi/sqrt(rowSums(xi^2))
    if (any(is.na(yi))){ # fill in columns with column mean when NA
      for(iii in 1:ncol(yi)){
        yi[is.na(yi[,iii]), iii] <- mean(yi[,iii], na.rm = TRUE)
      }
    }
    res <- NULL
    while( is.null(res) ) { # this will hang if there is a problem
      try(
        res <- kmeans(yi, centers)
      )
    } 
    #
    clus <- res$cluster
  }
  ## end clustering of whole data
  ## make judgement matrix
  rowNum <- nrow(t(allresults[[k]]$ordered_data))
  S <- matrix(0, rowNum, rowNum)
  for (j in 1:k) {
    X <- rep(0, rowNum)
    X[which(clus == j)] <- 1
    S <- S + X %*% t(X)
  }
  ## get perturbed probs
  rm <- allresults[[k]]$consensus_matrix
  rm <- as.matrix(rm)
  return(likelihood(S,rm))
}
