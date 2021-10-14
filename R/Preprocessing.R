
### ---------------------------- Method - 01 --------------------------------###
### Probabilistic Quotient Normalization (PQN).

PQN <- function(data) {
    #data<-log(data) #
        reference <- apply(data, 1, median)
        quotient <- data/reference
        quotient.median <- apply(quotient, 2, median)
        pqn.data <- t(t(data)/quotient.median)
        return(pqn.data)
}


### ---------------------------- Method - 02 --------------------------------###
### Cyclic Locally Weighted Regression (Cyclic Loess Normalization).

# load package unless it is already loaded
LOESS <- function(data) {
  loess.data <- affy::normalize.loess(data,
                                subset = 1:nrow(data),
                                epsilon = 10^-2,
                                maxit = 2,
                                log.it = FALSE, #
                                verbose = TRUE,
                                span = 0.75,
                                family.loess = "gaussian")
  return(loess.data)
}


### ---------------------------- Method - 03 --------------------------------###
### Contrast Normalization (Contrast).

# load package unless it is already loaded
CONTRAST <- function(data) {
  #---First adaption: Make the data matrix non-negative
  threshold = 1e-11
  data[data <= 0] <- threshold
  #---Apply normalization
  maffy.data <- affy::maffy.normalize(data,
                                subset = 1:nrow(data),
                                span = 0.75,
                                verbose = TRUE,
                                family = "gaussian",
                                log.it = FALSE) #
  #---Second adaption: Subtract 10% Quantile from each sample
  subtract <- function(x){
    t(t(x)-apply(x, 2, quantile, 0.1))
  }
  contrast.data <- subtract(maffy.data)

  rownames(contrast.data) <- rownames(data)
  colnames(contrast.data) <- colnames(data)

  return(contrast.data)
}


### ---------------------------- Method - 04 --------------------------------###
### Quantile Normalization (Quantile).

#load package unless it is already loaded
QUANTILE <- function(data) {
    #data<-log(data)##########################################################
    normalize.quantile <- get("normalize.quantiles",
                              envir = asNamespace("affy"))
  quantile.data <- normalize.quantile(data)
  rownames(quantile.data) <- rownames(data)
  colnames(quantile.data) <- colnames(data)
  return(quantile.data)
}


### ---------------------------- Method - 05 --------------------------------###
### Linear Baseline Normalization (Linear Baseline).

LINEAR <- function(data) {
    #data<-log(data)##########################################################
  linear.baseline <- apply(data,1,median)
  baseline.mean <- mean(linear.baseline)
  sample.means <- apply(data,2,mean)
  linear.scaling <- baseline.mean/sample.means
  linear.baseline.data <- t(t(data)*linear.scaling)
  return(linear.baseline.data)
}


### ---------------------------- Method - 06 --------------------------------###
### Non-Linear Baseline Normalization (Li-Wong Normalization).

#load package unless it is already loaded
LIWONG <- function(data) {
    #data<-log(data)##########################################################
  #---First step: Find baseline sample
  average.intensity <- apply(data, 2, mean)
  #R has an add way of rounding.
  median.number <- round(ncol(data)/2 + 0.1)
  #the additional 0.1 ensures that it rounds properly
  ordering <- order(average.intensity)
  median.sample.number <- ordering[median.number]
  median.sample <- data[, median.sample.number]
  #---Apply normalization
  liwong.data <- vector()
  for(i in 1:ncol(data)){
    liwong.model <- affy::normalize.invariantset(data = data[, i],
                                           ref = median.sample,
                                           prd.td = c(0.003,0.007))
    #the threshold of the rank-invariant set might need to be adjusted from case to case
    liwong.sample <- predict(liwong.model$n.curve$fit, data[,i])
    liwong.data <- cbind(liwong.data, liwong.sample$y)
  }
  return(liwong.data)
}

### ---------------------------- Method - 07 --------------------------------###
### Cubic Spline Normalization (Cubic Spline).
CUBIC <- function(data) {
  #data<-log(data)##########################################################
  #load package unless it is already loaded
  spline.data <- affy::normalize.qspline(data,
                                         samples = 0.1, # 0.02

                                         spline.method = "natural",

                                         target = apply(data, 1, mean))
  rownames(spline.data) <- rownames(data)
  colnames(spline.data) <- colnames(data)
  return(spline.data)
}
### ****************************************************************************
### code chunk number 03: Second group of methods - for between-variable variations.
### ****************************************************************************

### ---------------------------- Method - 08 --------------------------------###
### Auto Scaling (unit variance scaling).

AUTO <- function(data) {
  centered.data <- data - apply(data, 1, mean)
  scaling.auto <- apply(data, 1, sd)
  auto.data <- centered.data / scaling.auto
  return(auto.data)
}

# ----------- or using following manner.
#### Method - 01. Auto Scaling.
#auto.data <- DiffCorr::scalingMethods(obj, methods = "auto")
# ----------- End.

### ---------------------------- Method - 09 --------------------------------###
### Range Scaling.
RANGE <- function(data) {
  range.data <- DiffCorr::scalingMethods(data, methods = "range")
  return(range.data)
}

### ---------------------------- Method - 10 --------------------------------###
### Pareto Scaling.

PARETO <- function(data) {
  centered.data <- data - apply(data, 1, mean)
  scaling.pareto <- sqrt(apply(data, 1, sd))
  pareto.data <- centered.data / scaling.pareto
  return(pareto.data)
}

# ----------- or using following manner.
#### Method - 03. Pareto Scaling.
# pareto.data <- DiffCorr::scalingMethods(data, methods = "pareto")
# ----------- End.

### ---------------------------- Method - 11 --------------------------------###
### Vast Scaling.

VAST <- function(data) {
  vast.data <- DiffCorr::scalingMethods(data, methods = "vast")
  return(vast.data)
}


### ---------------------------- Method - 12 --------------------------------###
### Level Scaling.

LEVEL <- function(data) {
  level.data <- DiffCorr::scalingMethods(data, methods = "level")
  return(level.data)
}

### ---------------------------- Method - 13 --------------------------------###
### Variance Stabilization Normalization (VSN).

VSN <- function(data) {
  # load package unless it is already loaded
  vsn.model <- suppressMessages(vsn2(data))
  vsn.data <- suppressMessages(predict(vsn.model, data))
  return(vsn.data)
}


### ---------------------------- Method - 14 --------------------------------###
### Power Transformation.

POWER <- function(data) {
  power.data <- DiffCorr::scalingMethods(data, methods = "power")
  return(power.data)
}


### ---------------------------- Method - 15 --------------------------------###
### Log Transformation.

LOGTRAN <- function(data) {
  tmp_data <- apply(data, 1, function(x) log10(x) - mean(log10(x)))
  logtran.data <- data.frame(t(tmp_data), check.names = FALSE)
  return(logtran.data)
}


### ---------------------------- Method - 16 --------------------------------###
### MSTUS, the total signal 'MS total useful signal' (MSTUS).

MSTUS <- function(data) {
    data_sum <- matrix(colSums(data), nrow = 1)
    uni <- matrix(rep(1, nrow(data)), ncol = 1)
    area.uni <- uni %*% data_sum
    MSTUS <- data/area.uni
    return(MSTUS)
}

### ---------------------------- Method - 17 --------------------------------###
### Median. a method in R package metabolomics.

MEDIAN <- function(data) {
    #data <- log(data)
    inputdata<-data.frame(as.factor(rep("sample",ncol(data))),t(data))
    norm_med <- Normalise(inputdata, method = "median")
    median<-t(norm_med$output[,-1])
    return(median)
}


### ---------------------------- Method - 18 --------------------------------###
### Sum, a method in R package metabolomics.

SUM <- function(data) {
    #a <- log(data)
    a <- data
    inputdata<-data.frame(as.factor(rep("sample",ncol(a))),t(a))
    norm_sum <- Normalise(inputdata, method = "sum")
    sum<-t(norm_sum$output[,-1])
    return(sum)
}




### ---------------------------- Method - 19 --------------------------------###
### Mean, a method in R package metabolomics.

MEAN <- function(data) {
    #a <- log(data)
    a <- data
    inputdata<-data.frame(as.factor(rep("sample",ncol(a))),t(a))
    norm_mean <- Normalise(inputdata, method = "mean")
    mean<-t(norm_mean$output[,-1])
    return(mean)
}

### ---------------------------- Other methods --------------------------------###
###
### Data format:(1)Matrix: row is samples, column is metabolites.(The first column is the binary labels.)
###            (2)nc is the column order of QC metabolites or IS.


### ---------------------------- Method - 21 --------------------------------###
###  SIS method.
SIS <- function(data, nc) {
    # load package unless it is already loaded
    #library(metabolomics)
    norm_is <- Normalise(data, method = "is", refvec=data[, nc[1]])
    normdata_is <- norm_is$output[, -1]
    return(normdata_is)
}

### ---------------------------- Method - 22 --------------------------------###
###  NOMIS method.
NOMIS <- function(data, nc) {
    # load package unless it is already loaded
    #library(metabolomics)
    norm_nomis <- Normalise(data, method = "nomis", nc = nc)
    normdata_nomis <- norm_nomis$output[,-1]
    return(normdata_nomis)
}

### ---------------------------- Method - 23 --------------------------------###
###  CCMN method.
CCMN <- function(data, nc) {
    # load package unless it is already loaded
    #library(metabolomics)
    norm_ccmn <- Normalise(data, method = "ccmn", nc = nc, ncomp = 2)
    normdata_ccmn <- norm_ccmn$output[,-1]
    return(normdata_ccmn)
}

### ---------------------------- Method - 24 --------------------------------###
###  RUVRand method.
RUVRand <- function(data, nc) {
    # load package unless it is already loaded
    Y<-data.matrix(data[,-1])
    IS <- Y[,nc[1]]
    r<-numeric(dim(Y)[2])
    for(j in 1:length(r)){
        r[j]<-cor(IS,Y[,j])
    }
    ctl<-logical(length(r))
    ctl[which(r>round(quantile(r,0.7),2))]<-TRUE
    ruv <- MetNorm::NormalizeRUVRand(Y=Y,ctl=ctl,k=3,lambda=0.03, plotk = FALSE)
    #ruv <- NormalizeRUVRand(Y=Y,ctl=ctl,k=3, plotk = FALSE)
    normdata_ruvrand <- ruv$newY
    return(normdata_ruvrand)
}

### ---------------------------- END --------------------------------###

RUV2_dm <- function(data, nc, k) {
    # load package unless it is already loaded
    #library(metabolomics)
    treated.log <- data
    treated.group<-factor(treated.log[,1],levels=unique(treated.log[,1]))
    premat<-treated.log[which(treated.log[,1]=="case"),-1]
    postmat<-treated.log[which(treated.log[,1]=="control"),-1]
    ordFit<-LinearModelFit(datamat=data.matrix(postmat-premat),ruv2=TRUE, nc = nc, k = k,
                           factormat=matrix(1,nrow=nrow(postmat)), saveoutput=FALSE)
    return(ordFit)
}





### End.
### -------------------------------------------------------------------------###


### ----------------------------Cube Root Transformation --------------------------------###

Cube_root  <- function(data) {
  return(sign(data) * abs(data)^{1/3})
}





