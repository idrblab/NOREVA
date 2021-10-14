#' @title Time-course Metabolomic Study with dataset without QCSs and ISs.
#' @description this function enables the performance assessment of metabolomic data processing for time-course dataset
#' (without quality control sample and without internal standard) using four
#' independent criteria, and can comprehensively scan thousands of processing workflows and rank
#' all these workflows based on their performances (assessed from four different perspectives).
#' @param fileName Allows the user to indicate the NAME of peak table resulted from PrepareInuputFiles() (default = null).
#' @param SAalpha Allows the user to specify whether the input peak table satisfies the study assumption Alpha (SAalpha, all metabolites are assumed to be equally important) (default = “Y”).
#' “Y” denotes that the peak table satisfies the study assumption Alpha (SAalpha).
#' “N” denotes that the peak table does not satisfy the study assumption Alpha (SAalpha).
#' @param SAbeta Allows the user to specify whether the input peak table satisfies the study assumption Beta (SAbeta, the level of metabolite abundance is constant among all samples) (default = “Y”).
#' “Y” denotes that the peak table satisfies the study assumption Beta (SAbeta).
#' “N” denotes that the peak table does not satisfy the study assumption Beta (SAbeta).
#' @param SAgamma Allows the user to specify whether the input table satisfies study assumption Gamma (SAγ, the intensities of most metabolites are not changed under the studied conditions) (default = “Y”).
#' “Y” denotes that the peak table satisfies the study assumption Gamma (SAγ).
#' “N” denotes that the peak table does not satisfy the study assumption Gamma (SAγ).
#' @import DiffCorr affy vsn DT iterators
#' @import e1071 AUC impute MetNorm
#' @import ggsci timecourse multiROC dummies
#' @import ggplot2 ggord ggfortify usethis
#' @import ggrepel ggpubr sampling crmn
#' @rawNamespace import(limma, except=.__C__LargeDataObject)
#' @rawNamespace import(ropls, except=plot)
#' @importFrom grDevices dev.off png rainbow rgb colorRampPalette pdf
#' @importFrom graphics abline close.screen legend mtext par points screen split.screen symbols text title lines
#' @importFrom stats anova as.formula cor dnorm kmeans lm loess loess.control mad median model.matrix na.omit pf pnorm qnorm qt quantile rnorm runif sd var
#' @importFrom utils combn read.csv write.csv write.table
#' @usage nortimecoursenoall(fileName, SAalpha="Y", SAbeta="Y", SAgamma="Y")
#' @export nortimecoursenoall
#' @examples
#' library(NOREVA)
#' \donttest{timec_non_data <- PrepareInuputFiles(dataformat = 1,
#' rawdata = "Timecourse_without_QCSIS.csv")}
#' \donttest{nortimecoursenoall(fileName = timec_non_data,
#' SAalpha="Y", SAbeta="Y", SAgamma="Y")}


nortimecoursenoall <- function(fileName, SAalpha="Y", SAbeta="Y", SAgamma="Y"){

  cat("\n")
  cat("NOREVA is Running ...","\n")
  cat("\n")

  cat("*************************************************************************","\n")
  cat("Depending on the size of your input dataset","\n")
  cat("Several mintues or hours may be needed for this assessment","\n")
  cat("*************************************************************************","\n")
  cat("\n")

  cat("STEP 1: Prepare input file in standard formats of NOREVA", "\n")
  cat("\n")
  cat("STEP 2: The assumption(s) held as indicated by users","\n")
  cat("Study Assumption alpha: all proteins were equally important (Y/N): ", SAalpha, "\n")
  cat("Study Assumption beta: the level of protein abundance was constant among all samples (Y/N): ", SAbeta, "\n")
  cat("Study Assumption gamma: the intensity of the majority of proteins were unchanged (Y/N): ", SAgamma, "\n")

  cat("\n")
  cat("STEP 3: The criteira selected by users for this assessment","\n")
  cat("Criterion Ca: Reduction of Intragroup Variation", "\n")
  cat("Criterion Cb: Differential Metabolic Analysis", "\n")
  cat("Criterion Cc: Consistency in Marker Discovery", "\n")
  cat("Criterion Cd: Classification Accuracy", "\n")
  cat("\n")

  cat("NOREVA is Running ...","\n")
  cat("\n")
  #To ignore the warnings during usage
  options(show.error.messages=FALSE,echo=TRUE,keep.source.pkg=TRUE)
  #defaultW <- getOption("warn")
  options(warn=-1)

  consistency <-  function(fold = 3, top = 20) {
    folds <- fold

    DEG <- list()
    for (i in 1:folds) {

      set.seed(2)

      com.x <- t(test.fold[[i]][,-(1:2)])
      lab.ca <- as.factor(test.fold[[i]][,2])

      gnames <- rownames(com.x)

      ###different time points
      time.grp <- lab.ca

      ###times is the the number of time points
      times <- length(unique(lab.ca))

      ###A numeric vector or matrix corresponding to the sample sizes for all genes across different biological conditions, three classes in this case
      size <- rep(length(which(time.grp==unique(lab.ca)[1])), nrow(com.x))

      #out1 <- mb.long(fruitfly, times=12, reps=size, rep.grp = assay, time.grp = time.grp)
      out1 <- mb.long(com.x, times=times, reps=size, time.grp = time.grp)

      ### get marker ranking
      marker_ranking <- cbind(gnames, out1$HotellingT2)

      DEG[[i]] <- marker_ranking[order(as.numeric(marker_ranking[,2]), decreasing=T),1]


    }
    names(DEG) <- LETTERS[1:folds]

    top.n <-top # Extracting the top n genes.
    DEG.list <- DEG
    for (g in 1:length(DEG.list)) {
      DEG.list[[g]] <- DEG.list[[g]][1:top.n]
    }

    # Calculating consistency score:

    setlist <- DEG.list


    return(setlist) # consistense score

  }

  stabel.score <- function(repeats = 20, fold = 3, top = 10) {
    score <- 0
    for (r in 1:repeats) {
      score <- score + consistency(fold, top)
    }
    return(score/repeats)
  }
  ####----------------------------spssSkewKurtosis----------------------------------#######
  spssSkewKurtosis=function(x) {
    w=length(x)
    m1=mean(x)
    #print()
    m2=sum((x-m1)^2)
    m3=sum((x-m1)^3)
    m4=sum((x-m1)^4)
    s1=sd(x)

    skew = w*m3/(w-1)/(w-2)/s1^3

    sdskew=sqrt( 6*w*(w-1) / ((w-2)*(w+1)*(w+3)) )

    kurtosis=(w*(w+1)*m4 - 3*m2^2*(w-1)) / ((w-1)*(w-2)*(w-3)*s1^4)

    sdkurtosis=sqrt( 4*(w^2-1) * sdskew^2 / ((w-3)*(w+5)) )

    mat=matrix(c(skew, kurtosis, sdskew, sdkurtosis), 2,
               dimnames=list(c("skew","kurtosis"), c("estimate","se")))
    return(mat)
  }

  ###Imputation---------------------------------------------------------------------------------
  imput<-function(filter_data2,n){
     set.seed(3)
    matrix<-switch(
      n,
      t(ImputMean(filter_data2)),#
      t(ImputMedian(filter_data2)),#
      t(back(filter_data2)),#
      t(impute.knn(as.matrix(t(filter_data2)), k = 10,rng.seed = 1024)$data)#
    )
    return(matrix)
  }

  im_nam<-c(
    "MEI",
    "MDI",
    "HAM",
    "KNN"
  )

  ###Transformation---------------------------------------------------------------------------------
  trans<-function(data,n){
    matrix<-switch(
      n,
      Cube_root(data),
      log2(data),
      data
    )
    return(matrix)
  }

  t_nam<-c(
    "CUT",
    "LOG",
    "NON"
  )

  ###Normalization---------------------------------------------------------------------------------
  norm_sam<-function(train_data_t,n){
    matrix<-switch(
      n,
      train_data_t,#1
      PQN(train_data_t),#2
      LOESS(train_data_t),#3
      CONTRAST(train_data_t),#4
      QUANTILE(train_data_t),#5
      LINEAR(train_data_t),#6
      LIWONG(train_data_t),#7
      CUBIC(train_data_t),#8
      AUTO(train_data_t),#9
      RANGE(train_data_t),#10
      PARETO(train_data_t),#11
      VAST(train_data_t),#12
      LEVEL(train_data_t),#13
      VSN(train_data_t),#14
      POWER(train_data_t),#15
      MSTUS(train_data_t),#16
      SUM(train_data_t),#17
      MEDIAN(train_data_t),#18
      MEAN(train_data_t),#19
      EIGENMS(train_data_t, sampleLabel)#20
    )
    return(matrix)
  }

  n_nam<-c(
    "NON",
    "PQN",
    "LOE",
    "CON",
    "QUA",
    "LIN",
    "LIW",
    "CUB",
    "AUT",
    "RAN",
    "PAR",
    "VAS",
    "LEV",
    "VSN",
    "POW",
    "MST",
    "SUM",
    "MED",
    "MEA",
    "EIG"
  )

  norm_no <- "NON"

  #normalization for with IS---------------------------------------------------------------------------------
  #-----------------------------------------------------------------
  ###################################################Step-2 READ DATA
  data_q <- fileName
  #data_q <- read.csv("Data_Time-course_without_QCS_and_ISs.csv", header = TRUE,stringsAsFactors = FALSE)
  data_q<-data_q[,-2]
  data2 <- data_q[, -c(1:2)]
  data4 <- data_q[, 1:2]
  data2[data2 == 0] <- NA # the zero value has been replaced by NA.

  #filtering-----------------------------------------------------------------------------
  col_f <- apply(data2, 2, function(x) length(which(is.na(x)))/length(x))
  if (length(which(col_f >0.2))==0){
    data2_f <- data2
  }else {
    data2_f <- data2[, -which(col_f > 0.2)]
  }
  filter_data2 <- data2_f

  #---------------------------------------------------------------------------------
  aftetable<-list()
  normal_data <- list()
  train_data_metabolite3 <- NULL

  #############################################################
  lengthA = SAalpha
  lengthB = SAbeta
  lengthC = SAgamma

  normal1 <- NA
  normal3 <- NA
  normal5 <- NA
  ### N1 --- All methods ---###
  if(any(lengthA == "Y") && any(lengthB == "Y") && any(lengthC == "Y")){
    normal1 <- c(2:8,16:20)
    normal2 <- c(1,9:13,15)
    normal3 <- c(9:13,15)
    normal4 <- c(1,2:8,16:20)
    normal5 <- c(1,14)
  }

  ### N2 --- Assumption A: N; Assumption A: N; Assumption A: N ---###
  if(any(lengthA == "N") && any(lengthB == "N") && any(lengthC == "N")){
    normal5 <- c(1,14,20)
  }

  ### N3 --- Assumption A: Y; Assumption A: N; Assumption A: N ---###
  if(any(lengthA == "Y") && any(lengthB == "N") && any(lengthC == "N")){
    normal1 <- c(9,13,11,10,12,15)
    normal2 <- 1
  }

  ### N4 --- Assumption A: N; Assumption A: Y; Assumption A: N ---###
  if(any(lengthA == "N") && any(lengthB == "Y") && any(lengthC == "N")){
    normal1 <- c(4,8,3,6,7,2,5)
    normal2 <- 1
  }

  ### N5 --- Assumption A: N; Assumption A: N; Assumption A: Y ---###
  if(any(lengthA == "N") && any(lengthB == "N") && any(lengthC == "Y")){
    normal1 <- c(19,18,17,16)
    normal2 <- 1
  }

  ### N6 --- Assumption A: Y; Assumption A: Y; Assumption A: N ---###
  if(any(lengthA == "Y") && any(lengthB == "Y") && any(lengthC == "N")){
    normal1 <- c(2:8)
    normal2 <- c(9:13,15)

    normal3 <- c(9:13,15)
    normal4 <- c(2:8)
  }

  ### N7 --- Assumption A: Y; Assumption A: N; Assumption A: Y ---###
  if(any(lengthA == "Y") && any(lengthB == "Y") && any(lengthC == "N")){
    normal1 <- c(16:19)
    normal2 <- c(9:13,15)

    normal3 <- c(9:13,15)
    normal4 <- c(16:19)
  }

  ### N8 --- Assumption A: N; Assumption A: Y; Assumption A: Y ---###
  if(any(lengthA == "N") && any(lengthB == "Y") && any(lengthC == "Y")){
    normal1 <- c(2:8)
    normal2 <- c(16:19)

    normal3 <- c(16:19)
    normal4 <- c(2:8)
  }

  #################################################################
  sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
  #for (i in as.numeric(impt)){
  k.l <- foreach (i = 1:4,.combine = "c") %do%{
    tryCatch({
    afterimpute.table <- NULL
    imput_m <- imput(filter_data2,i)

    imputed_data <- cbind(data4, imput_m)
    afterimpute.table <- imputed_data
    getLabel <-as.character(data_q[, 2])
    data1 <- afterimpute.table
    train_data_t <- t(data1[, -(1:2)])

    dd <- spssSkewKurtosis(train_data_t)
    dd1 <- as.data.frame(dd)
    DD2 <- dd1$estimate/dd1$se
    mat <- t(matrix(DD2, 2, dimnames=list(c("skew_zscore","kurtosis_zscore"))))
    mat1 <- as.data.frame(mat)
    szscore <- mat1$skew_zscore
    kzscore <- mat1$kurtosis_zscore

    dd2 <- as.data.frame(t(matrix(dd1$estimate,2, dimnames=list(c("skew1","kurtosis1")))))
    skew <- dd2$skew1
    kurtosis <- dd2$kurtosis1
    if ((any(nrow(data1) < 50)&&any(abs(szscore) <= 1.96)) || (any(nrow(data1) < 50)&&any(abs(kzscore) <= 1.96))){
      #cat("No transformation")
      tform <- 3
    }else if((any(nrow(data1) <300)&&any(nrow(data1) > 50)&&any(abs(szscore) <= 3.29)) || (any(nrow(data1) <300)&&any(nrow(data1) > 50)&&any(abs(kzscore) <= 3.29))){
      #cat("No transformation")
      tform <- 3
    }else if((any(nrow(data1) >300)&&any(abs(skew) <= 2)) || (any(nrow(data1) >300)&&any(abs(kurtosis) <= 7))){
      #cat("No transformation")
      tform <- 3
    }else{
      #cat("Use transformation")
      tform <- c(1,2)
    }
      }, error=function(e){})
    #for (j in as.numeric(trsf)){
    foreach (j = tform) %do%{

      train_data_Transformation3<-try(trans(train_data_t,j))
      if(class(train_data_Transformation3)=="try-error")
      return(NA)
      #train_data_t <- train_data_Transformation3
      train_data_Transformation3[is.infinite(data.matrix(train_data_Transformation3))]<-NA
      sampleLabel <- as.character(afterimpute.table[, 2])

      imputed_data <- cbind(data4, t(train_data_Transformation3))
      after.table <- imputed_data

      return(list(i,j,train_data_Transformation3=train_data_Transformation3,after.table=after.table,sampleLabel=sampleLabel,afterimpute.table=afterimpute.table))
        }}
      ####
      ## calculate
      cluster <- makeCluster(parallel::detectCores()-1, type = "SOCK") ;cluster %>% registerDoSNOW ; time = proc.time()

      k.norm12 <-  foreach(n=k.l,.packages=c("foreach","dplyr","vsn") ,.combine = "c") %dopar% {
        i <- n[[1]]
        #q <- n[[2]]
        j <- n[[2]]
        train_data_Transformation3 <- n$train_data_Transformation3
        afterqc.table <- n$afterqc.table
        sampleLabel <- n[[5]]
        afterimpute.table<-n$afterimpute.table

        k.n1 <- list()

      if(any(!is.na(normal1))){
        k.n1 <- foreach(k = normal1,.combine = "c") %do% {
          set.seed(3)
      #for ( k in normal1){
        train_data_Preprocess3 <-try(norm_sam(train_data_Transformation3,k))
        if(class(train_data_Preprocess3)=="try-error")
          return(NA)

        foreach (h =normal2)%do%{
        #for(h in normal2){

          train_data_metabolite3<-try(norm_sam(train_data_Preprocess3,h))
          if(class(train_data_metabolite3)=="try-error")
            return(NA)
          set.seed(3)
          normalized_data3 <- try(t(train_data_metabolite3))
          if(class(normalized_data3)=="try-error")
            return(NA)

          eva.data3 <- cbind(afterimpute.table[, 1:2], normalized_data3)
          eva.data3 <- eva.data3[, -1]
          eva.data3 <- as.data.frame(eva.data3)
          rownames(eva.data3) <- afterimpute.table[, 1]
          colnames(eva.data3)[1] <- "Group"
          eva_data3<-eva.data3
          k.a <- list()
          k.b <- list()
          k.a[[paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")]] <- eva_data3

          k.b[[paste(im_nam[i],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- after.table
          return(list(k.a,k.b))

          #normal_data[[paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")]] <- eva_data3
          #eva_data3 <- NULL
          #aftetable[[paste(im_nam[i],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- after.table
          #save(normal_data,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
        }
      }}
        k.n2 <- list()
      if(any(!is.na(normal3))){
      #for(k in normal3){
        k.n2 <-  foreach(k = normal3,.combine = "c")%do%{
        set.seed(3)
        train_data_Preprocess3 <-try(norm_sam(train_data_Transformation3,k))
        if(class(train_data_Preprocess3)=="try-error")
          return(NA)
        foreach(h = normal4)%do%{
          set.seed(3)
          train_data_metabolite3<-try(norm_sam(train_data_Preprocess3,h))
          if(class(train_data_metabolite3)=="try-error")
            return(NA)

          normalized_data3 <- try(t(train_data_metabolite3))
          if(class(normalized_data3)=="try-error")
            return(NA)

          eva.data3 <- cbind(afterimpute.table[, 1:2], normalized_data3)
          eva.data3 <- eva.data3[, -1]
          eva.data3 <- as.data.frame(eva.data3)
          rownames(eva.data3) <- afterimpute.table[, 1]
          colnames(eva.data3)[1] <- "Group"
          eva_data3<-eva.data3

          k.a <- list()
          k.b <- list()
          k.a[[paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")]] <- eva_data3
          k.b[[paste(im_nam[i],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- after.table
          return(list(k.a,k.b))

          #normal_data[[paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")]] <- eva_data3

          #eva_data3 <- NULL

          #aftetable[[paste(im_nam[i],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- after.table
          #save(normal_data,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
        }
      }}
      k.n3 <- list()
      if(any(!is.na(normal5))){
        k.n3 <- foreach(k = normal5) %do%{
          set.seed(3)

        train_data_metabolite3<-try(norm_sam(train_data_Transformation3,k))
        if(class(train_data_metabolite3)=="try-error")
          return(NA)

        normalized_data3 <- try(t(train_data_metabolite3))
        if(class(normalized_data3)=="try-error")
          return(NA)

        eva.data3 <- cbind(afterimpute.table[, 1:2], normalized_data3)
        eva.data3 <- eva.data3[, -1]
        eva.data3 <- as.data.frame(eva.data3)
        rownames(eva.data3) <- afterimpute.table[, 1]
        colnames(eva.data3)[1] <- "Group"
        eva_data3<-eva.data3

        k.a <- list()
        k.b <- list()
        k.a[[paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k]),"]", sep=""), sep="+")]] <- eva_data3
        k.b[[paste(im_nam[i],t_nam[j],n_nam[k],sep="+")]] <- after.table
        return(list(k.a,k.b))
        #normal_data[[paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k]),"]", sep=""), sep="+")]] <- eva_data3

        #eva_data3 <- NULL
        #aftetable[[paste(im_nam[i],t_nam[j],n_nam[k],sep="+")]] <- after.table
        #save(normal_data,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
      }}
      return(c(k.n1,k.n2,k.n3))
    }
print(proc.time()-time)
stopCluster(cluster)
  sink()

  Fpmad<-list()
  Fpurity<-list()
  Fscore<-list()
  Fauc<-list()
  mmm <-144
  ########### all list name ###########
  newname=c()
  for (i in 1:4){
    for (j in tform){
      if(any(!is.na(normal1))){
        for ( k in normal1){
          for(h in normal2){
            namebind <- paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")
            newname=rbind(newname,namebind)
          }
        }}
      if(any(!is.na(normal3))){
        for(k in normal3){
          for(h in normal4){
            namebind <- paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")
            newname=rbind(newname,namebind)
          }
        }}
      if(any(!is.na(normal5))){
        for(k in normal5){
          namebind <- paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_nam[k]),"]", sep=""), sep="+")
          newname=rbind(newname,namebind)
        }}
    }
  }

  k.test <- k.norm12 %>%sapply(.,function(x)  x[[1]] ,USE.NAMES = TRUE)
  kk <- k.norm12[!is.na(k.test)]
  normal_data  <- kk  %>%sapply(.,function(x) x[[1]] , USE.NAMES = TRUE)
  aftetable   <- kk  %>%sapply(.,function(x)  x[[2]] ,USE.NAMES = TRUE)
  sink()

  # save(normal_data,after.table,aftetable,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
  #save(normal_data,aftetable,newname,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
  save(kk,newname,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")

  #save(normal_data,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
  #load("./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
  options(show.error.messages=FALSE,echo=TRUE,keep.source.pkg=TRUE)
  #defaultW <- getOption("warn")
  options(warn=-1)

  dir.create(paste0("OUTPUT-NOREVA-Criteria.Ca"))
  dir.create(paste0("OUTPUT-NOREVA-Criteria.Cb"))
  dir.create(paste0("OUTPUT-NOREVA-Criteria.Cc"))
  dir.create(paste0("OUTPUT-NOREVA-Criteria.Cd"))

  #################################################################
  #options(show.error.messages=FALSE,echo=FALSE,keep.source.pkg=TRUE)
  #defaultW <- getOption("warn")
  #options(warn=-1)
  #nanmes_right<-names(normal_data)
  # clean useless data
  rm(afterqc.table,afterqc_table,data_q,data1,data2,
     data2_f,data2_QC,fileName,filter_data2,imput_m2,
     imput_m2_t,imputed_data,k.l,k.list,k.norm12,
     k.test,aftetable,normal_data,samFile,samPeno,sampleData,
     sampleData_rev,sampleLabel,sampleList,train_data_t,train_data_Transformation)

  ################################Step 2
  opts <- list(progress=function(n) setTxtProgressBar(txtProgressBar(min=1, max=length(kk), style=3), n))
  cluster <- makeCluster(parallel::detectCores()-1, type = "SOCK") ;cluster%>% registerDoSNOW ; time = proc.time() #
  # k.test <- foreach::foreach (i = 1:length(normal_data),.packages=c("reshape")) %dopar% k.pp(i,afterqc.table,normal_data,nanmes_right,newname) # length(normal_data)
  k.test <- foreach::foreach (k.input = iterators::iter(kk),.options.snow=opts,.packages=c("reshape","ggplot2"),.combine = "rbind") %dopar% {  # ,.errorhandling = c("pass")
    k.step2.name <- names(k.input[[1]])
    k.normal_data <- k.input[[1]][[1]]
    k.aftetable <- k.input[[2]][[1]]
    # ============================================准备一些函数=========================================####
    consistency_M <- function(fold = 3, top = 20) {

      folds <- fold

      DEG <- list()
      for (i in 1:folds) {

        com.x <- test_data[test.fold[[i]], ]

        com.x[sapply(com.x, simplify = 'matrix', is.infinite)] <- 0

        X_matrix <- as.data.frame(com.x[,-1])
        y_label <- as.factor(com.x[,1])

        pos_filter <- OPLSDA_C(X_matrix, y_label)

        DEG[[i]] <- pos_filter
      }
      names(DEG) <- LETTERS[1:folds]

      top.n <- top # Extracting the top n genes.
      DEG.list <- DEG
      for (g in 1:length(DEG.list)) {
        DEG.list[[g]] <- DEG.list[[g]][1:top.n]
      }

      # Calculating consistency score:

      setlist <- DEG.list


      return(setlist) # consistense score

    }
    spssSkewKurtosis=function(x) {
      w=length(x)
      m1=mean(x)
      m2=sum((x-m1)^2)
      m3=sum((x-m1)^3)
      m4=sum((x-m1)^4)
      s1=sd(x)
      skew = w*m3/(w-1)/(w-2)/s1^3
      sdskew=sqrt( 6*w*(w-1) / ((w-2)*(w+1)*(w+3)) )
      kurtosis=(w*(w+1)*m4 - 3*m2^2*(w-1)) / ((w-1)*(w-2)*(w-3)*s1^4)
      sdkurtosis=sqrt( 4*(w^2-1) * sdskew^2 / ((w-3)*(w+5)) )
      mat=matrix(c(skew, kurtosis, sdskew, sdkurtosis), 2,
                 dimnames=list(c("skew","kurtosis"), c("estimate","se")))
      return(mat)
    }


    #name <- nanmes_right
    ####
    id <- which(newname==k.step2.name, arr.ind = TRUE)
    id <- as.data.frame(id)
    id1 <- id$row

  #for (mmm in 1:length(normal_data)){
    #mmm=1
    #name <- nanmes_right

    ####
    #id <- which(newname==k.step2.name, arr.ind = TRUE)
    #id <- as.data.frame(id)
    #id1 <- id$row
    #####

    ###Fpmad-------------------------------------------------------------------------
    n_data <- as.data.frame(k.normal_data,col.names=NULL)
    n_data <- as.matrix(n_data)
    if(sum(is.na(n_data))<length(n_data)/3){
      eva_data3<-as.data.frame(k.normal_data,col.names=NULL)
    }else{return(NA)}

    pmad3N.log <- eva_data3
    pmad3N <- try(PMAD(pmad3N.log))
    if(class(pmad3N)=="try-error")
    { return(NA) }

   # Fpmad[names(k.normal_data)]<-mean(pmad3N)

    names(k.aftetable)[1] <- "Group" # input afterqc.table
    ##########################################################################################
    #### input
    ##########################################################################################
    pmad3R.log <- k.aftetable[,-1]
    pmad3R <- try(PMAD(pmad3R.log))
    if(class(pmad3R)=="try-error")
    { return(NA) }

    C3 <- cbind(pmad3R, pmad3N); colnames(C3) <- c("Before", "After")

    cat(paste("Assessing Method" , paste(id1,"/",length(newname),":",sep=""), k.step2.name),"\n")

    cat("   Criterion Ca (reduction of intragroup variation) ...","\n")
    #dir.create(paste0("OUTPUT-NOREVA-Criteria.Ca"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Ca/Criteria.Ca-", k.step2.name,".pdf",sep=""))
    #library(reshape2)
    C3new <- as.data.frame(C3)
    melt1C3 <- cbind(C3new, "name" = rownames(C3new))
    melt2C3 <- melt(melt1C3, id.vars = "name")
    colnames(melt2C3) <- c("name", "beforeafter", "value")
    try(print(ggplot(melt2C3, aes(x = melt2C3$beforeafter, y = melt2C3$value, color = melt2C3$beforeafter)) +
                geom_violin(width=0.5,size=1.5) +
                scale_color_manual(values=c("#fbbc05","#800080")) +
                geom_boxplot(color=c("#fbbc05","#800080"), size=1.5, width=0.1)+
                theme_bw() +
                theme(panel.grid.major =element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.title.y=element_blank(),
                      axis.title.x=element_blank(),
                      legend.position="none"
                )
    ))
    dev.off()
    ###2. Fpurity-------------------------------------------------------------------------

    data_kmeans <- eva_data3
    data_kmeans[sapply(data_kmeans, simplify = 'matrix', is.na)] <- 0

    eva_data3_f<-data_kmeans
    del_col <- NULL
    for (m in 1:dim(eva_data3_f)[2]){
      if (sum(eva_data3_f[,m]==0)==nrow(eva_data3_f)) {
        del_col <- c(del_col,m)
      }
    }

    if (is.null(del_col)){
      eva_data3_f <- eva_data3_f
    }else{
      eva_data3_f <- eva_data3_f[,-del_col]
    }

    com.x <- t(eva_data3_f[,-1])
    lab.ca <- as.factor(eva_data3_f[,1])
    gnames <- rownames(com.x)

    ###different time points
    time.grp <- lab.ca

    ###times is the the number of time points
    times <- length(unique(lab.ca))

    ###A numeric vector or matrix corresponding to the sample sizes for all genes across different biological conditions, three classes in this case
    size <- rep(length(which(time.grp==unique(lab.ca)[1])), nrow(com.x))

    sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
    out1 <- try(mb.long(com.x, times=times, reps=size, time.grp = time.grp))
    if(class(out1)=="try-error")
    { return(NA) }

    marker_ranking <- cbind(gnames, out1$HotellingT2)

    DEG <- marker_ranking[order(as.numeric(marker_ranking[,2]), decreasing=T),]

    #data_kmeans <- as.data.frame(eva_data3_f[, c("Group",DEG[(1:50),1])])
    data_kmeans <- as.data.frame(eva_data3_f[, c("Group",DEG[(1:20),1])])

    clusters <- length(unique(data_kmeans[, 1]))
    obj_kmeans <- try(kmeans(data_kmeans[,-1], centers = clusters, nstart = 1,iter.max = 3))
    if(class(obj_kmeans)=="try-error")
    { return(NA) }

    groups <- factor(data_kmeans[, 1], levels = unique(data_kmeans[,1]))
    unique.groups <- levels(groups)
    cols <- 1:length(unique(data_kmeans[, 1]))
    box_cols <- NULL

    for (ii in 1:length(data_kmeans[, 1])) {

      box_cols[ii] <- cols[match(data_kmeans[, 1][ii],unique.groups)]
    }

    true_label<-box_cols

    pre<-obj_kmeans$cluster
    tru<-true_label
    tmatrix<-table(pre,tru)

    label<-tru
    result<-pre

    accuracy<-try(purity(result,label))
    if(class(accuracy)=="try-error")
    { return(NA) }

   # Fpurity[names(k.normal_data)]<-accuracy

    unique.groups <- levels(as.factor(data_kmeans[,1]))
    T_number <- length(unique.groups)
    cols <- rainbow(length(unique.groups))

    box_cols <- c(rep(NA, length(rownames(data_kmeans))))

    for (ii in 1:length(data_kmeans[, 1])) {
      box_cols[ii] <- cols[match(data_kmeans[, 1][ii],unique.groups)]
    }

    data_kmeans <- as.data.frame(data_kmeans)
    data_kmeans$Color<-box_cols

    data_kmeans$T_label <- data_kmeans[, 1]

    sink()

    cat("   Criterion Cb (differential metabolic analysis) ...","\n")

    sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
    #dir.create(paste0("OUTPUT-NOREVA-Criteria.Cb"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cb/Criteria.Cb-",k.step2.name,".pdf",sep=""))
    classcolor <- NA
    for (i in 1:length(unique(obj_kmeans$cluster))){
      majority <- obj_kmeans$cluster[obj_kmeans$cluster==i]
      majority1 <- names(majority)
      majority2 <- data_kmeans[majority1,"Color"]
      majority3 <- unique(majority2)[which.max(tabulate(match(majority2, unique(majority2))))]
      majority4 <- as.character(majority3)
      classcolor[i] <- majority4
    }

    try(print(autoplot(obj_kmeans,
                       data = data_kmeans,
                       frame = TRUE,shape=1)+
                geom_text_repel(aes(label=data_kmeans$T_label),color=data_kmeans$Color,size=4,family="serif")+
                theme(legend.position="none",panel.background = element_blank(),axis.title.x=element_text(angle=0, size=12,color="black"),axis.text.x=element_text(angle=0, size=13,color="black"),axis.title.y=element_text(size=12,color="black"),axis.text.y=element_text(size=13,color="black"),panel.border = element_rect(fill='transparent', color='black'))+
                geom_point(color=data_kmeans$Color,size=3)+
                #scale_fill_manual(values = c("#BC4D70", "#00B1C1", "#55A51C")) +
                #scale_color_manual(values = c("#BC4D70", "#00B1C1", "#55A51C"))
                scale_fill_manual(values = classcolor) +
                scale_color_manual(values = classcolor)))
    dev.off()

    ###3. Consistency -------------------------------------------------------- #
    test_data <- eva_data3
    test_data[sapply(test_data, simplify = 'matrix', is.na)] <- 0

    Sample<-rownames(test_data)
    test_data<-cbind(Sample,test_data)

    for(iii in 1:dim(test_data)[1]){
      test_data$Sample_label[iii]<-strsplit(as.character(test_data$Sample),"T")[[iii]][1]
    }
    test_data[1:5,1:5]
    Sample_labels<-unique(test_data$Sample_label)

    filter_label1 <- sample(Sample_labels,round(length(Sample_labels)/3),replace = FALSE)
    filter_label2<-sample(Sample_labels[-match(filter_label1,Sample_labels)],round(length(Sample_labels)/3),replace = FALSE)
    filter_label3<-Sample_labels[-c(match(filter_label1,Sample_labels), match(filter_label2,Sample_labels))]

    test.fold <- list()
    group1<-test_data[test_data$Sample_label %in% filter_label1,]
    group1$Sample_label<-NULL
    test.fold[[1]] <- group1

    group2<-test_data[test_data$Sample_label %in% filter_label2,]
    group2$Sample_label<-NULL
    test.fold[[2]] <- group2

    group3<-test_data[test_data$Sample_label %in% filter_label3,]
    group3$Sample_label<-NULL
    test.fold[[3]] <- group3

    # DEG.list<-try(consistency(3, 50))
    # DEG.list<-try(consistency(3, 20))
    tryCatch({
    DEG.list<- consistency(3, 10)}, error=function(e){})
    #if(class(DEG.list)=="try-error")
    #{ next }

    CW_value <- try(CWvalue(DEG.list,Y=(ncol(eva_data3)-1),n=length(DEG.list[[1]])))
    if(class(CW_value)=="try-error")
    { return(NA) }

    # Fscore[names(k.normal_data)]<-CW_value

    sink()

    cat("   Criterion Cc (consistency in marker discovery) ...","\n")

    sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))

    #dir.create(paste0("OUTPUT-NOREVA-Criteria.Cc"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cc/Criteria.Cc-",k.step2.name,".pdf",sep=""))
    try(print(plot(eulerr::venn(DEG.list),fills = list(fill = c("white", "white","white")),
                   labels = list(col = "black", font = 2),
                   edges = list(col = c("#800080", "#4285f4", "#fbbc05"), lwd=4),
                   quantities = TRUE)))
    dev.off()
    ### -- 4. AUC value --------------------------------------------------------------------------------- #
    DEG <- marker_ranking[order(as.numeric(marker_ranking[,2]), decreasing=T),]

    ######################################################################

    set.seed(3)

    # NB ROC PLOT will change for each new random noise component (jitter)
    X_matrix <- eva_data3[,-1]
    y_label <- as.factor(eva_data3[,1])
    X_matrix[sapply(X_matrix, simplify = 'matrix', is.na)] <- 0

    data_multiROC <- cbind(y_label, X_matrix)
    #data_multiROC <- cbind(data_multiROC[,-1], paste("Label_", data_multiROC[, 1], sep = ""))
    data_multiROC <- cbind(data_multiROC[,-1], paste("Label_", data_multiROC[, 1], sep = ""), stringsAsFactors=TRUE)
    colnames(data_multiROC)[ncol(data_multiROC)] <- "Label"

    X_matrix <- data_multiROC
    y_label <- data_multiROC[,"Label"]

    #x <- as.data.frame(X_matrix[, c(DEG[(1:50),1], "Label")])
    x <- as.data.frame(X_matrix[, c(DEG[(1:20),1], "Label")])
    y <- y_label
    y<- as.factor(x[, length(x)])
    folds <- 5
    test.fold <- split(sample(1:length(y)), 1:folds) #ignore warning

    for (mm in 1:5) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]
      dim(test_df)
      svmmodel <- try(e1071::svm(Label ~ ., data = train_df, type = "C-classification", probability = TRUE))
      if(class(svmmodel)=="try-error")
      { return(NA) }

      svm_pred <- try(predict(svmmodel, test_df,probability = TRUE))
      if(class(svm_pred)=="try-error")
      { return(NA) }

      svm_pred <- data.frame(attr(svm_pred, "probabilities"))
      colnames(svm_pred) <- paste(colnames(svm_pred), "_pred_SVM")

      true_label <- dummies::dummy(test_df$Label, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste(colnames(true_label), "_true")

      final_df_2 <- cbind(true_label, svm_pred)
      final_df2_t<-t(final_df_2)

      final_df2_t_la<-rownames(final_df2_t)
      final_df2_t_c<-as.data.frame(cbind(final_df2_t_la,final_df2_t))

      if (mm==1){final_df_m<-final_df2_t_c}
      else{
        final_df_m<-merge(final_df_m,final_df2_t_c,by="final_df2_t_la",all =T)
      }

    }
    final_df_m[is.na(final_df_m)]<-0
    rownames(final_df_m)<-final_df_m[,1]

    final_df_m_t<-as.data.frame(t(final_df_m[,-1]))
    roc_res <- try(multiROC::multi_roc(final_df_m_t, force_diag = T))
    if(class(roc_res)=="try-error")
    { return(NA) }

    plot_roc_df <- try(multiROC::plot_roc_data(roc_res))
    if(class(plot_roc_df)=="try-error")
    { return(NA) }

    plot_roc_df_mic<-subset(plot_roc_df,plot_roc_df$Group=="Micro")

    AUC_mic<-round(unique(plot_roc_df_mic$AUC),3)


  #  Fauc[names(k.normal_data)]<-AUC_mic


    sink()

    cat("   Criterion Cd (classification accuracy) ...","\n")
    cat("\n")

    #dir.create(paste0("OUTPUT-NOREVA-Criteria.Cd"))

    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cd/Criteria.Cd-",k.step2.name,".pdf",sep=""))
    plot(x=c(0,1),y=c(0,1),col="lightgrey",pch=16,bg="yellow",type = 'l',xlim=c(0,1),ylim=c(0,1),lwd=1,xlab="1-Specificity",ylab="Sensitivity")
    lines(x=1-plot_roc_df_mic$Specificity,plot_roc_df_mic$Sensitivity,col="red",pch=16,bg="yellow",xlim=c(0,1),ylim=c(0,1),lwd=2,xlab="WEEK",ylab="STUDE")

    dev.off()
    return(c(k.step2.name,mean(pmad3N),accuracy,CW_value,AUC_mic))

  }

  print(proc.time()-time)
  stopCluster(cluster)
  save(k.test,file="./step2_data.Rdata")

  k.result1 <- k.test%>%.[,-1]%>%apply(.,2,as.numeric)%>% data.table(.,id=k.test%>%.[,1]) %>% .[is.na(V1)==F]
  k.result2 <- k.result1[,.(V1,V2,V3,V4)] %>% as.data.frame()
  colnames(k.result2) <- c("Precision","Cluster_accuracy","Reproducibility","Classification")
  rownames(k.result2) <- k.result1$id
  result2 <- k.result2

  #需要加载的包data.table,tidyverse,plyr,doSNOW,foreach,parallel

  # =========================================排名csv及排名热图的输出======================================####

 # result<-dplyr::bind_rows("Precision"=unlist(Fpmad),"Cluster_accuracy"=unlist(Fpurity),"Reproducibility"=unlist(Fscore),"Classification"=unlist(Fauc),.id = "id")
 # result1<-t(result)
 # colnames(result1)<-result1["id",]
 # result2<-result1[-1,]
 # result2<-data.frame(result2,check.names=FALSE)

  for(i in 1:dim(result2)[2]){result2[,i]=as.numeric(as.character(result2[,i]))}

  Rank<-apply(result2, 2, function(x){rank(-x,ties.method="min",na.last = "keep")})

  if(length(grep("Precision",colnames(Rank)))==1){
    Rank[,"Precision"]<-rank(as.numeric(as.character(result2[,"Precision"])),ties.method="min",na.last = "keep")
  }else{
    Rank<-Rank
  }

  Rank_revision<-apply(Rank, 2, function(x){x[is.na(x)]<-nrow(Rank);return(x)})
  Ranksum0<-apply(Rank_revision, 1, sum)
  Rankres0<-cbind("OverallRank"=Ranksum0,Rank_revision)

  zuihou0<-cbind("Rank"=Rankres0,"Value"=result2)

  zuihou1<-zuihou0[order(Rankres0[,"OverallRank"],decreasing = FALSE),]
  zuihou1[,1]<-rank(zuihou1[,1],ties.method="min")
  zuihou2<-zuihou1
  zuihou3 <- zuihou2
  zuihou3 <- round(zuihou3,4)
  colnames(zuihou3) <- c("Overall-Rank","Criteria.Ca-Rank","Criteria.Cb-Rank","Criteria.Cc-Rank","Criteria.Cd-Rank","Criteria.Ca-Value","Criteria.Cb-Value","Criteria.Cc-Value","Criteria.Cd-Value")

  ##########picture######################################
  data_color<-as.data.frame(zuihou2[,-c(1:5)])

  data_color["Value.Precision"][data_color["Value.Precision"]>=0.7]<-1
  data_color["Value.Precision"][data_color["Value.Precision"]<0.7&data_color["Value.Precision"]>=0.3]<-8
  data_color["Value.Precision"][data_color["Value.Precision"]<0.3]<-10

  data_color["Value.Cluster_accuracy"][data_color["Value.Cluster_accuracy"]>=0.8]<-10
  data_color["Value.Cluster_accuracy"][data_color["Value.Cluster_accuracy"]<0.8&data_color["Value.Cluster_accuracy"]>=0.5]<-8
  data_color["Value.Cluster_accuracy"][data_color["Value.Cluster_accuracy"]<0.5]<-1

  data_color["Value.Reproducibility"][data_color["Value.Reproducibility"]>=0.3]<-10
  data_color["Value.Reproducibility"][data_color["Value.Reproducibility"]<0.3&data_color["Value.Reproducibility"]>=0.15]<-8
  data_color["Value.Reproducibility"][data_color["Value.Reproducibility"]<0.15]<-1

  data_color["Value.Classification"][data_color["Value.Classification"]>=0.9]<-10
  data_color["Value.Classification"][data_color["Value.Classification"]<0.9&data_color["Value.Classification"]>=0.7]<-8
  data_color["Value.Classification"][data_color["Value.Classification"]<0.7]<-1

  data_color_m<-as.data.frame(data_color)
  Ranksum_color<-apply(data_color_m, 1, sum)
  data_color_m01<-cbind( "rank_color"=Ranksum_color,data_color_m)
  data_color_m02<-data_color_m01[order(data_color_m01[,"rank_color"],decreasing =T),]
  data_color_m<-data_color_m02

  row<-rownames(data_color_m)
  nfina<-nchar(row[1])
  nstart<-nchar(row[1])-6
  result <- substring(row, nstart,nfina)

  data_heat<-data_color_m[,-1]
  colnames(data_heat) <- c("Criterion Ca: Reduction of Intragroup Variation","Criterion Cb: Differential Metabolic Analysis","Criterion Cc: Consistency in Marker Discovery","Criterion Cd: Classification Accuracy")

  rank_result <- zuihou3[match(row.names(data_heat),row.names(zuihou3)),]
  rank_result[,1] <- 1:nrow(rank_result)

  write.csv(rank_result,file = "./OUTPUT-NOREVA-Overall.Ranking.Data.csv")

  #options(warn = defaultW)
  cat("\n")
  cat("*************************************************************************","\n")
  cat("Congratulations! Assessment Successfully Completed!","\n")
  cat("Thanks for Using NOREVA. Wish to See You Soon ...","\n")
  cat("*************************************************************************","\n")
  cat("\n")
  #return(rank_result)
}
