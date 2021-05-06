loess.fastlo <- function (y, yhat, degree=NULL, weights=NULL, control=NULL) {

  tfit <- loess(y - yhat ~ yhat, degree=degree, weights=weights, control=control)
  smooth <- predict(tfit, yhat[2:( length(yhat) - 1 )])

  return(smooth)
}
##============================ fastlo ========================#
fastlo <- function(x, subset, maxit=2, mfun=NULL, log.it = FALSE, verbose=TRUE,
                   epsilon = 0.01, MM=F, parallel=FALSE, ...) {

  x[is.na(data.matrix(x))]<-0

  #####
  # "quick help"
  #####
  if ( is.null(x) )
  {
    msg <- paste(
      "\n##########",
      "\n# fastlo()  - Version 1.1.0 - Summer 2004 (estimated)",
      "\n#   x       = matrix of values (rows are 'probes' and columns are 'chips')",
      "\n#   subset  = index into a subset of rows to use in the normalization (default is 'all rows')",
      "\n#   maxit   = maximum number of iterations (default is 3)",
      "\n#   mfun    = function to use for estimating yhat (default is 'yhat <- rowMeans(y - smooth)')",
      "\n#   log.it  = do we want to log the PM values before fitting them (default is TRUE)",
      "\n#   verbose = how much information to report back during process (default is TRUE)",
      "\n#   epsilon = convergence criteria (default is 0.01)",
      "\n#   parallel = boolean indicating whether algorithm should be run in parallel using SNOW package (default is FALSE)",
      "\n##########\n",
      sep="")
    #cat(msg)
    invisible(return(FALSE))
  }

  lc<-loess.control(surface="interpolate",statistics="approximate",
                    trace.hat="approximate",cell=0.2,iterations=4)

  #### HAD TO ADD FOR R
  if (class(x) == 'Plob') {
    if (MM) y<-rbind(x@pm,x@mm)
    else y<-x@pm
  }

  else {
    if (class(x) == 'matrix') y<-x
    else { #cat("x must be a matrix or Probe Level Object","\n" )
      return() }
  }

  dtemp <- dimnames(y)     #save the dimnames of y for later
  dimnames(y) <- NULL      # this speeds things up

  nchip <- ncol(y)
  nspot <- nrow(y)

  if(!missing(subset)){
    if(length(subset)==1) {   #if a subset exists take a sample
      subset <- sample(1:nspot, subset)
    }
    else {
      if(any(subset < 1 | subset > nspot))
        stop("Invalid rows in the subset argument")
    }
  }
  else {    ##when subset = missing takes all of the rows
    subset <- 1:nspot
  }

  if(log.it) y <- logb(y,2)

  # initialize matrices with same dimensions as of y but filled with zeros.
  smooth <- matrix(0., ncol=nchip, nrow=nspot)
  old.smooth <- 0

  w <- c(0,rep(1,length(subset)),0) ##weights for loess fit
  cl <- NULL

    cl <- parallel::makeCluster(nchip)


  for (iter in 1:maxit) {
    #if(verbose) cat("Iteration", iter, "Chip")

    # TODO: parallelize rowMeans as well?
    if (is.null(mfun)) yhat <- rowMeans(y - smooth)  # this is faster
    else               yhat <- apply(y - smooth, 1, mfun)

    # To avoid NA's due to extrapolation, the lowess has to use
    #   include the largest and smallest yhat
    temp <- order(yhat)
    index <- c(temp[1], subset, temp[nspot])

    #####
    # Here's where we're going to try to go parallel.
    # We need to wrap the loess model code in a new function
    # that we can pass to parCapply().  Then we can do each
    # loess fit on a separate node in the cluster.
    #####
      smooth <- matrix(data=parallel::parCapply(cl, y[index,], loess.fastlo, yhat=yhat[index], degree=1, weights=w, control=lc), ncol=nchip)
      for(j in 1:nchip) {
        # separate function loess.fastlo allows parallelization
        smooth[,j] <- loess.fastlo(y=y[index,j], yhat=yhat[index], degree=1, weights=w, control=lc)
      }


    change <- max(colMeans((smooth[subset,] - old.smooth)^2))
    old.smooth <- smooth[subset,]
    #if (verbose) cat("\n   Finished, change = ", format(change),"\n")
    if (change <= epsilon) break
  }

  if(!is.null(cl)) {
    parallel::stopCluster(cl)
  }
  #if (verbose) cat("\n")

  ynorm <- y - smooth    # normalized y
  dimnames(ynorm) <- dtemp
  if(log.it) ynorm<-2^ynorm
  if(MM) {
    pm<-ynorm[1:(nspot/2), ]
    mm<-ynorm[((nspot/2)+1):nspot, ]
    return(pm,mm) }
  else
    return(ynorm)
}
