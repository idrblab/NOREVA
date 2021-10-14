#' @title shiftCor
#' @description shiftCor provides the QC-RLS correction for large scale metabolomics. See the details at the following References.
#' @param samPeno a file with  the meta information.
#' @param samFile a file with  the expression information.
#' @param Frule The cut-off value for missing value filter function.
#' @param QCspan The smoothing parameter which controls the bias-variance tradeoff. if the QCspan is set at '0', the generalised cross-validation will be performed to avoid overfitting the observed data.
#' @param degree Lets you specify local constant regression (i.e., the Nadaraya-Watson estimator, degree=0), local linear regression (degree=1), or local polynomial fits (degree=2, the default).
#' @param imputeM The parameter for imputation method.(i.e., nearest neighbor averaging, "KNN"; minimum values for imputed variables, "min", median values for imputed variables (Group dependent) "median").
#' @examples
#' library(statTarget)
#' datpath <- system.file("extdata",package = "statTarget")
#' samPeno <- paste(datpath,"MTBLS79_sampleList.csv", sep="/")
#' samFile <- paste(datpath,"MTBLS79.csv", sep="/")
#' shiftCor(samPeno,samFile,Frule = 0.8,QCspan = 0.75, degree = 2,imputeM = "KNN")
#' @references statTarget: a streamlined tool for signal drift correction
#' and interpretations of quantitative mass spectrometry-based omics data.
#' Luan H, Ji F, Chen Y, Cai Z. 2018, Analytica Chimica Acta.
#' @export
#' @import statTarget
shiftCor <- function(samPeno,samFile,Frule = 0.8,QCspan = 0.75, degree = 2,imputeM = "KNN"){
  #imputeM = "median","KNN" or "min"
  #degree = "0","1","2"
  #samPeno <- read.csv(samPeno, header=TRUE, check.names = F,stringsAsFactors = F)
  samPeno <- as.data.frame(samPeno)
  #samFile <- read.csv(samFile,header=FALSE, check.names = F,stringsAsFactors = F)
  samFile <- t(samFile)
  #colnames(samFile) <- samFile[1,]
  samFile <- as.data.frame(samFile[-1,])
  #rownames(samFile) <- samFile$name
  samFP <- samFile

  #if(sum(is.na(samPeno$class)) <=0 ){
    #message(date(),"There were not QC sample in your data!")
    #return(para)
  #}

  #################LOESS data NA########
  samFP <- as.matrix(samFP)
  samFP[samFP<=0] <- NA
  #cat(date(), "\nstatTarget: shiftCor start... \nEvaluation of missing value...")
  #message(date(),"The number of NA value in Data Profile before QC-RLSC: ",
          #sum(is.na(samFP) | as.matrix(samFP) <= 0))
  imsamFP <- samFP
  #############Filter miss value###################

  FilterMV = function(m,degre) {
    dx <- c()
    for(i in 1:ncol(m)){
      freq <- as.vector(tapply(m[,i], degre, function(x){sum(is.na(x) | as.matrix(x) <= 0)/length(x)}))
      if(sum(freq > Frule) > 0) dx <- c(dx , i)
    }
    if(length(dx) >0) m <- m[,-dx]
    return(m)
  }
  classF <- as.factor(samPeno$class)
  classF = addNA(classF)
  imsamFPF = FilterMV(imsamFP,classF)
  #Frule_warning= paste("The number of vaiables including", Frule*100, "% of missing value :",sep = " ")
  #message(date(),Frule_warning," ", dim(imsamFP)[2]-dim(imsamFPF)[2])
  imsamFP = imsamFPF
  ##############impute missing value#################
  #cat(date(), "\nImputation start...\n")
  if(imputeM == "KNN"){
    #require(impute)
    mvd <- impute::impute.knn(imsamFP[,1:ncol(imsamFP)])
    inputedData <- mvd$data
  }else if(imputeM == "min"){
    inputedData <- apply(imsamFP[,2:ncol(imsamFP)],2,function(y){
      y[is.na(y) | y<=0] <- min(y[y>0],na.rm = TRUE)
      y})
    #inputedData <- t(inputedData)
  }else if(imputeM == "median"){
    missvalue <- function(x,group) {
      x[is.na(x) == TRUE ] <- 0
      group = as.factor(as.numeric(group))
      for (i in 1:dim(x)[1]){
        for(j in 2:dim(x)[2]){
          if(x[i,j] == 0 | sum(is.na(x[i,j])) >= 1){
            #x[i,j][is.na(x[i,j]) == TRUE ] <- 0
            x[i,j] <- tapply(as.numeric(x[,j]),group,median)[group[i]]
          }
        }
      }
      return(x)
    }
    inputedData = missvalue(imsamFP,classF)
    inputedData = inputedData[,-1]
  }
  #message(date(), "The number of NA value in Data Profile after the initial imputation: ",
          #sum(is.na(inputedData) | as.matrix(inputedData) <= 0))

  if(sum(is.na(inputedData) | as.matrix(inputedData) <= 0) > 0)
  {
    inputedData[inputedData<=0] <- NA
    mvd2 <- impute::impute.knn(inputedData[,1:ncol(inputedData)])
    inputedData <- mvd2$data
    #message(date(), "The number of NA value in Data Profile after the second imputation (KNN): ",
            #sum(is.na(inputedData) | as.matrix(inputedData) <= 0))
  }

  #cat(date(), "\nImputation Finished!\n")


  dat <- as.matrix(t(inputedData))
  numX <- 1:dim(dat)[2]

  if(QCspan > 0){
    loessFit=function(x,y,QCspan,degree){
      cn <- colnames(x)
      ####Check########
      st_QC<-grep("QC",cn[1])
      ed_QC<-grep("QC",cn[length(cn)])
      if(length(st_QC)==0)
      {
        stop("Wrong: the first sample must be QC sample; please check ......");
      }
      if(length(ed_QC)==0)
      {
        stop("Wrong: the sample at the end of sequence must be QC sample; please check ......");
      }
      qcid <- grep("QC",cn)
      for(i in 1:dim(x)[1]){
        loe <- stats::loess(x[i,qcid]~qcid,span=QCspan,degree = degree)
        yf <- stats::predict(loe,y)
        x[i,] <- as.numeric(x[i,])/yf
      }
      loessDat = x
    }
    loessDat <- loessFit(x = dat,y=numX,QCspan = QCspan,degree = degree)
  }else if(QCspan <= 0){
    #message(date(),"\nWarning: The QCspan was set at '0'. \nThe LOESS based generalised cross-validation was used to avoid overfitting the observed data")
    autoFit <- function(xl,y){
      cn <- colnames(xl)
      ####Check########
      st_QC<-grep("QC",cn[1])
      ed_QC<-grep("QC",cn[length(cn)])
      if(length(st_QC)==0)
      {
        stop("Wrong: the first sample must be QC sample; please check ......");
      }
      if(length(ed_QC)==0)
      {
        stop("Wrong: the sample at the end of sequence must be QC sample; please check ......");
      }
      qcid <- grep("QC",cn)
      for(i in 1:dim(xl)[1]){
        loe1 <- loess(xl[i,qcid]~qcid)
        loe2 <- loe1
        sp <- c(seq(0.05,0.75,0.01))
        CVspan <-c()
        for(j in 1:length(sp)){
          mod <- stats::update(loe1, span = sp[j])
          CVspan[j] = loessGCV(mod)[["gcv"]]
        }
        minG <- as.matrix(data.frame(sp,CVspan))
        minG[!is.finite(minG)] <- max(minG[,2],na.rm = TRUE)
        minspan <- minG[which.min(minG[,2]),1]
        minspan
        #spanNew <- spanCV(loe)$minimum
        loeN <- stats::update(loe2, span = minspan)
        yf <- predict(loeN,y)
        xl[i,] <- as.numeric(xl[i,])/yf
      }
      loessDat = xl
      #return(loessDat)
    }
    loessDat <- autoFit(xl = dat,y = numX)
  }
  ###############


  if(sum(is.na(loessDat) | as.matrix(loessDat) <= 0) > 0)
  {
    loessDat[loessDat<=0] <- NA
    mvd2 <- impute::impute.knn(loessDat[,1:ncol(loessDat)])
    loessDat <- mvd2$data
    #message(date(),"The number of NA value in Data Profile after Loess Correction (KNN): ",
            #sum(is.na(loessDat) | as.matrix(loessDat) <= 0))
  }
  #dirout.uni = paste(getwd(), "/statTarget/", sep = "")
  #dir.create(dirout.uni)
  #dirout.w = paste(getwd(), "/statTarget/shiftCor", sep="")
  #dir.create(dirout.w)
  #dirout.Bs = paste(getwd(), "/statTarget/shiftCor/Before_shiftCor", sep="")
  #dir.create(dirout.Bs)
  #dirout.As = paste(getwd(), "/statTarget/shiftCor/After_shiftCor", sep="")
  #dir.create(dirout.As)
  ########### Out plot of each peak ############
  #for(i in 1 :dim(dat)[1]){
  #  loplot(dat,loessDat,i)
  #}
  ###############Raw output###########

  raw_temp <- cbind(samPeno,inputedData)
  nam_qc <- rownames(raw_temp)
  QC_temp_raw <- grep("QC",nam_qc)
  QC_temp_raw <- raw_temp[c(QC_temp_raw),]
  raw_temp_qc <- QC_temp_raw[,-c(3,4)]
  rownames(raw_temp_qc) <- NULL
  #RSD30_CV=paste("shift_QC_raw",".csv", sep="")
  #write.csv(raw_temp_qc, paste(dirout.Bs, RSD30_CV, sep="/"))
  #cat(date(), "\nCalculation of CV distribution of raw peaks (QC)...\n")
  #dirout.rs = paste(getwd(), "/statTarget/shiftCor/RSDresult/RSDdist_raw", sep="")
  #dir.create(dirout.rs)
  #setwd(dirout.rs)
  #RsdCal(paste(dirout.Bs, RSD30_CV, sep="/"))
  ###############Loess output###########

  lo_temp <- cbind(samPeno,t(loessDat))
  nam_qc <- rownames(lo_temp)
  QC_temp <- grep("QC",nam_qc)
  QC_temp <- lo_temp[-c(QC_temp),]
  lo_temp_sam <- QC_temp[,-c(2,4)]
  rownames(lo_temp_sam) <- NULL
  #RSD30_CV=paste("shift_sample_loess",".csv", sep="")
  #write.csv(lo_temp_sam, paste(dirout.As, RSD30_CV, sep="/"))
  return(lo_temp_sam)

  QC_temp <- grep("QC",nam_qc)
  QC_temp <- lo_temp[c(QC_temp),]
  lo_temp_qc <- QC_temp[,-c(3,4)]
  rownames(lo_temp_qc) <- NULL
  #RSD30_CV=paste("shift_QC_loess",".csv", sep="")
  #write.csv(lo_temp_qc, paste(dirout.As, RSD30_CV, sep="/"))
  #cat(date(), "\nCalculation of CV distribution of corrected peaks (QC)...\n")
  #dirout.rs = paste(getwd(), "/statTarget/shiftCor/RSDresult/RSDdist_cor", sep="")
  #dir.create(dirout.rs)
  #setwd(dirout.rs)
  #RsdCal(paste(dirout.As, RSD30_CV, sep="/"))

  #QC_temp <- grep("QC",nam_qc)
  #QC_temp <- lo_temp[c(QC_temp),]
  lo_temp_all <- lo_temp[,-c(2,4)]
  rownames(lo_temp_all) <- NULL
  #RSD30_CV=paste("shift_all_loess",".csv", sep="")
  #write.csv(lo_temp_all, paste(dirout.As, RSD30_CV, sep="/"))
  #write.table(lo_temp,"loess_out.csv",sep="/")
  #write.table(lo_temp,"",sep="\t",quote=F)
  #cat(date(), "\nCorrection Finished!\n")
  #cat(date(),"\nAll done!\n")
  ##################Loess Plot########################
}
