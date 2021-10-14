#' @title Multi-class (N>1) Metabolomic Study with dataset with Quality Control Samples (QCSs) and the corresponding data of golden standards for performance evaluation using Criterion e.
#' @description this function enables the performance assessment of metabolomic data processing
#' for multi-class dataset (with quality control sample but without internal standard)
#' using five independent criteria, and can comprehensively scan thousands of processing
#' workflows and rank all these workflows based on their performances (assessed from four different perspectives).
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
#' @param GS Allows the user to indicate the name of the file that contains the spike-in compounds (default = null).
#' The file should be in a .csv format, which provides the concentrations of spike-in compounds.
#' @import DiffCorr affy vsn DT  doSNOW
#' @import e1071 AUC impute MetNorm parallel
#' @import ggsci timecourse multiROC dummies
#' @import ggplot2 ggord ggfortify usethis
#' @import ggrepel ggpubr sampling crmn foreach dplyr
#' @rawNamespace import(limma, except=.__C__LargeDataObject)
#' @rawNamespace import(ropls, except=plot)
#' @importFrom grDevices dev.off png rainbow rgb colorRampPalette pdf
#' @importFrom graphics abline close.screen legend mtext par points screen split.screen symbols text title
#' @importFrom stats anova as.formula cor dnorm kmeans lm loess loess.control mad median model.matrix na.omit pf pnorm qnorm qt quantile rnorm runif sd var
#' @importFrom utils combn read.csv write.csv write.table
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @usage normulticlassqcallgs(fileName, GS, SAalpha="Y", SAbeta="Y", SAgamma="Y")
#' @export normulticlassqcallgs
#' @examples
#' library(NOREVA)
#' \donttest{multi_qc_data <- PrepareInuputFiles(dataformat = 1,
#' rawdata = "Multiclass_with_QCS.csv")}
#' \donttest{normulticlassqcallgs(fileName = multi_qc_data,
#' GS = "Multiclass_with_QCS_GoldenStandard.csv", SAalpha="Y", SAbeta="Y", SAgamma="Y")}

normulticlassqcallgs <- function(fileName, GS, SAalpha="Y", SAbeta="Y", SAgamma="Y"){

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
  cat("Criterion Ce: Level of correspondence between normalized and reference data", "\n")
  cat("\n")

  cat("NOREVA is Running ...","\n")
  cat("\n")

  #To ignore the warnings during usage
  options(show.error.messages=FALSE,echo=TRUE,keep.source.pkg=TRUE)
  #defaultW <- getOption("warn")
  options(warn=-1)

  #imputation---------------------------------------------------------------------------------
  consistency_M <-  function(fold = 3, top = 20) {

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

    top.n <-top # Extracting the top n genes.
    DEG.list <- DEG
    for (g in 1:length(DEG.list)) {
      DEG.list[[g]] <- DEG.list[[g]][1:top.n]
    }

    # Calculating consistency score:

    setlist <- DEG.list


    return(setlist) # consistense score

  }

  imput<-function(filter_data2,n){
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

  #Transformation---------------------------------------------------------------------------------
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

  #normalization---------------------------------------------------------------------------------
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
  ##----------------------------spssSkewKurtosis----------------------------------#######
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

  q_nam<-c(
    "NWE",
    "LLR",
    "LPF"
  )
  #-----------------------------------------------------------------

  ###################################################Step-2 调用数据
  data_q <- fileName
  sampleList <- as.data.frame(data_q[, 1:4])
  names(sampleList) <- c("sample", "batch", "class", "order")
  sampleData <- as.data.frame(data_q[, -(1:4)])

  sampleLabel_00 <- as.character(data_q[, 1])

  data2 <- sampleData
  data2_QC <- data_q[which(is.na(data_q[,3])),-(1:4)]

  #filtering-----------------------------------------------------------------------------
  col_f <- apply(data2_QC, 2, function(x) length(which(is.na(x)))/length(x))

  if (length(which(col_f > 0.2))==0){
    data2_f <- data2
  }else {
    data2_f <- data2[, -which(col_f >0.2)]
  }

  col_r <- apply(data2_QC, 2, function(x) sd(x)/mean(x))

  if (length(which(col_r > 0.3))==0){
    filter_data2 <- data2_f
  }else {
    filter_data2 <- data2_f[, -which(col_r > 0.3)]
  }
  dim(filter_data2)
  #---------------------------------------------------------------------------------

  normal_data <- list()
  aftetable<-list()
  eva_data3 <- NULL
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

  #sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
# define function
  f.tform <- function(train_data_t) {
    dd <- spssSkewKurtosis(train_data_t)
    dd1 <- as.data.frame(dd)
    DD2 <- dd1$estimate / dd1$se
    mat <-
      t(matrix(DD2, 2, dimnames = list(c(
        "skew_zscore", "kurtosis_zscore"
      ))))
    mat1 <- as.data.frame(mat)
    szscore <- mat1$skew_zscore
    kzscore <- mat1$kurtosis_zscore

    dd2 <-
      as.data.frame(t(matrix(dd1$estimate, 2, dimnames = list(
        c("skew1", "kurtosis1")
      ))))
    skew <- dd2$skew1
    kurtosis <- dd2$kurtosis1
    if ((any(nrow(data1) < 50) &&
         any(abs(szscore) <= 1.96)) ||
        (any(nrow(data1) < 50) && any(abs(kzscore) <= 1.96))) {
      #cat("No transformation")
      tform <- 3
    } else if ((any(nrow(data1) < 300) &&
                any(nrow(data1) > 50) &&
                any(abs(szscore) <= 3.29)) ||
               (any(nrow(data1) < 300) &&
                any(nrow(data1) > 50) && any(abs(kzscore) <= 3.29))) {
      #cat("No transformation")
      tform <- 3
    } else if ((any(nrow(data1) > 300) &&
                any(abs(skew) <= 2)) ||
               (any(nrow(data1) > 300) && any(abs(kurtosis) <= 7))) {
      #cat("No transformation")
      tform <- 3
    } else{
      #cat("Use transformation")
      tform <- c(1, 2)
    }
    return(tform)
  }

  f.eva <-  function(afterqc.table, x) {
    eva.data3 <- cbind(afterqc.table[, 1:2], normalized_data3)
    eva.data3 <- eva.data3[,-1]
    eva.data3 <- as.data.frame(eva.data3)
    rownames(eva.data3) <- afterqc.table[, 1]
    colnames(eva.data3)[1] <- "Group"
    return(eva.data3)
  }
  # calculate parameter
  k.l <-    foreach (i = 1:4,.combine = "c") %do%{

    afterqc.table<-NULL
    imput_m2 <- try(imput(filter_data2,i))
    if(class(imput_m2)=="try-error") return(NA)


    # running the QC processing.
    imput_m2_t<-t(imput_m2)
    sampleData_rev <- as.data.frame(cbind(rownames(imput_m2_t), imput_m2_t))
    colnames(sampleData_rev)<- c("name",sampleLabel_00)
    samPeno <- sampleList
    samFile <- sampleData_rev

    foreach (q = 1:3,.combine = "c") %do%{
      tryCatch({
        degree <- q-1
        afterqc_table <- shiftCor(samPeno, samFile, Frule = 0.8, QCspan = 0.75, degree = degree)

        afterqc.table<-afterqc_table
        getLabel <-as.character(afterqc.table[, 2])
        data1 <- afterqc.table
        train_data_t <- t(data1[, -(1:2)])

        tform <- f.tform(train_data_t) # 绮剧畝浜嗕竴涓猣unction
        #colnames(train_data_t) <- NULL
      }, error=function(e){})
      #for (j in as.numeric(trsf)){
      foreach (j = tform) %do%{
        tryCatch({
          train_data_Transformation <- trans(train_data_t,j)

          sampleLabel <- as.character(afterqc.table[, 2])

          imputed_data <- cbind(afterqc_table[,1:2],t(train_data_Transformation))

          afterqc_table<-imputed_data
          afterqc.table<-afterqc_table
        }, error=function(e){})

        return(list(i,q,j,train_data_Transformation=train_data_Transformation,afterqc.table=afterqc.table,sampleLabel=sampleLabel))
      }
    }
  }


  ## calculate

  cluster <- makeCluster(parallel::detectCores()-1, type = "SOCK") ;cluster %>% registerDoSNOW ; time = proc.time()

  k.norm12 <-  foreach(n=k.l,.packages=c("foreach","dplyr","vsn") ,.combine = "c") %dopar% {
    i <- n[[1]]
    q <- n[[2]]
    j <- n[[3]]
    train_data_Transformation <- n$train_data_Transformation
    afterqc.table <- n$afterqc.table
    sampleLabel <- n[[6]]

    k.n1 <- list()
    if(any(!is.na(normal1))){
      k.n1 <- foreach(k = normal1,.combine = "c") %do% {
        tryCatch({
          train_data_Preprocess3 <- norm_sam(train_data_Transformation,k)

        }, error=function(e){})

        foreach (h =normal2)%do%{

          tryCatch({
            normalized_data3 <- norm_sam(train_data_Preprocess3,h) %>% t
            eva_data3<-f.eva(afterqc.table,normalized_data3)
          }, error=function(e){})

          k.a <- list()
          k.b <- list()
          k.a[[paste(im_nam[i], q_nam[q], t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")]] <- eva_data3
          k.b[[paste(im_nam[i],q_nam[q],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- afterqc.table

          return(list(k.a,k.b))
        }
      }}

    k.n2 <- list()
    if(any(!is.na(normal3))){
      k.n2 <-  foreach(k = normal3,.combine = "c")%do%{

        tryCatch({
          train_data_Preprocess3 <- norm_sam(train_data_Transformation,k)
        }, error=function(e){})

        foreach(h = normal4)%do%{
          train_data_metabolite3<-try(norm_sam(train_data_Preprocess3,h))
          if(class(train_data_metabolite3)=="try-error")
          { h <- h+1 }

          normalized_data3 <- try(t(train_data_metabolite3))
          if(class(normalized_data3)=="try-error") return(NA)

          eva_data3<-f.eva(afterqc.table,normalized_data3)

          k.a <- list()
          k.b <- list()
          k.a[[paste(im_nam[i], q_nam[q], t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")]] <- eva_data3
          k.b[[paste(im_nam[i],q_nam[q],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- afterqc.table
          return(list(k.a,k.b))
        }
      }}

    k.n3 <- list()
    if(any(!is.na(normal5))){
      k.n3 <- foreach(k = normal5) %do%{
        train_data_metabolite3<-try(norm_sam(train_data_Transformation,k))
        if(class(train_data_metabolite3)=="try-error") return(NA)

        normalized_data3 <- try(t(train_data_metabolite3))
        if(class(normalized_data3)=="try-error") return(NA)

        eva_data3<-f.eva(afterqc.table,normalized_data3)

        k.a <- list()
        k.b <- list()
        k.a[[paste(im_nam[i], q_nam[q], t_nam[j], paste("[",paste(n_nam[k]),"]", sep=""), sep="+")]] <- eva_data3
        k.b[[paste(im_nam[i],q_nam[q],t_nam[j],n_nam[k],sep="+")]] <- afterqc.table
        return(list(k.a,k.b))
      }}

  #sink()


    return(c(k.n1,k.n2,k.n3))
  }
  print(proc.time()-time)
  stopCluster(cluster)

  # generation full set mathod
  k.list <- k.l %>% sapply(function(x) c(x[[1]],x[[2]],x[[3]]) ) %>% t %>% data.table
  newname <- NULL
  for (i in 1:4){
    for (q in 1:3){
      jk <- k.list %>% .[V1==i&V2==q,V3]
      for (j in jk){
        if(any(!is.na(normal1))){
          for ( k in normal1){
            for(h in normal2){
              namebind <- paste(im_nam[i], q_nam[q], t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")
              newname=rbind(newname,namebind)
            }
          }}
        if(any(!is.na(normal3))){
          for(k in normal3){
            for(h in normal4){
              namebind <- paste(im_nam[i], q_nam[q], t_nam[j], paste("[",paste(n_nam[k],n_nam[h],sep="-"),"]", sep=""), sep="+")
              newname=rbind(newname,namebind)
            }
          }}
        if(any(!is.na(normal5))){
          for(k in normal5){
            namebind <- paste(im_nam[i], q_nam[q], t_nam[j], paste("[",paste(n_nam[k]),"]", sep=""), sep="+")
            newname=rbind(newname,namebind)
          }}
      }
    }}
  # transform output to single var
  k.test <- k.norm12 %>%sapply(.,function(x)  x[[1]] ,USE.NAMES = TRUE)
  kk <- k.norm12[!is.na(k.test)]
  normal_data  <- kk  %>%sapply(.,function(x) x[[1]] , USE.NAMES = TRUE)
  aftetable   <- kk  %>%sapply(.,function(x)  x[[2]] ,USE.NAMES = TRUE)
  #sink()
  #print(proc.time()-time)
  save(normal_data,aftetable,newname,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
  #print(proc.time()-time)

  nanmes_right<-names(normal_data)

  # load("./OUTPUT-NOREVA-All.Normalized.Data.Rdata")

  # Fpmad<-list()
  # Fpurity<-list()
  # Fscore<-list()
  # Fauc<-list()

  options(show.error.messages=FALSE,echo=TRUE,keep.source.pkg=TRUE)
  #defaultW <- getOption("warn")
  options(warn=-1)

  dir.create(paste0("OUTPUT-NOREVA-Criteria.Ca"))
  dir.create(paste0("OUTPUT-NOREVA-Criteria.Cb"))
  dir.create(paste0("OUTPUT-NOREVA-Criteria.Cc"))
  dir.create(paste0("OUTPUT-NOREVA-Criteria.Cd"))
  dir.create(paste0("OUTPUT-NOREVA-Criteria.Ce"))
  # clean useless data
  # rm(afterqc.table,afterqc_table,data_q,data1,data2,
  #    data2_f,data2_QC,fileName,filter_data2,imput_m2,
  #    imput_m2_t,imputed_data,k.l,k.list,k.norm12,
  #    k.test,kk,samFile,samPeno,sampleData,
  #    sampleData_rev,sampleLabel,sampleList,train_data_t,train_data_Transformation)

  nanmes_right<-names(normal_data)

  ################################Step 2
  opts <- list(progress=function(n) setTxtProgressBar(txtProgressBar(min=1, max=length(normal_data), style=3), n))
  cluster <- makeCluster(parallel::detectCores()-1, type = "SOCK") ;cluster%>% registerDoSNOW ; time = proc.time() #
  k.test <- foreach::foreach (mmm = 1:length(normal_data), .options.snow=opts,.packages=c("reshape","ggplot2"),.combine = "rbind") %dopar% { # ,.errorhandling = c("pass")
    # ============================================准备一些函数=========================================####
    consistency_M <-  function(fold = 3, top = 20) {

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

      top.n <-top # Extracting the top n genes.
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


    name <- nanmes_right
    ####
    id <- which(newname==name[mmm], arr.ind = TRUE)
    id <- as.data.frame(id)
    id1 <- id$row
    #####################################################################################################
    ###1、Fpmad-------------------------------------------------------------------------
    #####################################################################################################

    ##########################################################################################
    #### input
    ##########################################################################################
    n_data <- as.data.frame(normal_data[mmm],col.names=NULL) # input-normal_data
    n_data <- as.matrix(n_data)
    if(sum(is.na(n_data))<length(n_data)/3){
      eva_data3<-as.data.frame(normal_data[mmm],col.names=NULL)

    }else{return(NA)}

    pmad3N.log <- eva_data3
    pmad3N <- try(PMAD(pmad3N.log))
    if(class(pmad3N)=="try-error")
    { return(NA) }

    ##########################################################################################
    #################### output
    ##########################################################################################
    # Fpmad[names(normal_data[mmm])]<-mean(pmad3N) # output
    ##########################################################################################
    #################### output
    ##########################################################################################
    # names(afterqc.table)[1] <- "Group" # input afterqc.table
    names(aftetable[[mmm]])[1] <- "Group" # input afterqc.table
    ##########################################################################################
    #### input
    ##########################################################################################
    pmad3R.log <- aftetable[[mmm]][,-1]
    pmad3R <- try(PMAD(pmad3R.log))
    if(class(pmad3R)=="try-error")
    { return(NA) }

    C3 <- cbind(pmad3R, pmad3N); colnames(C3) <- c("Before", "After")

    cat(paste("Assessing Method" , paste(id1,"/",length(newname),":",sep=""), name[mmm]),"\n")

    cat("   Criterion Ca (reduction of intragroup variation) ...","\n")

    dir.create(paste0("OUTPUT-NOREVA-Criteria.Ca"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Ca/Criteria.Ca-",names(normal_data[mmm]),".pdf",sep=""))
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
    #####################################################################################################
    ###2. Fpurity-------------------------------------------------------------------------
    #####################################################################################################
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

    X_matrix<-eva_data3_f[,-1]
    Group <- as.factor(eva_data3_f$Group)

    #sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
    tryCatch({
      pos_filter <- OPLSDA_test(X_matrix, Group, cutoff = 0.8)}, error=function(e){})
    ##sink()

    DEG <- cbind(Group,X_matrix[, pos_filter])
    data_kmeans<-DEG

    clusters <- length(unique(data_kmeans[, 1]))
    obj_kmeans <- try(kmeans(data_kmeans[,-1], centers = clusters, nstart = 1))
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
    ##########################################################################################
    #################### output
    ##########################################################################################
    # Fpurity[names(normal_data[mmm])]<-accuracy

    unique.groups <- levels(as.factor(data_kmeans[,1]))
    T_number <- length(unique.groups)

    cols <- rainbow(length(unique.groups))
    #cols <- colorRampPalette(c("#55a51c", "#ea7125", "#8f2bbc","#00b1c1"))(length(unique.groups))

    box_cols <- c(rep(NA, length(rownames(data_kmeans))))

    for (ii in 1:length(data_kmeans[, 1])) {
      box_cols[ii] <- cols[match(data_kmeans[, 1][ii],unique.groups)]
    }

    data_kmeans <- as.data.frame(data_kmeans)
    data_kmeans$Color<-box_cols

    data_kmeans$T_label <- data_kmeans[, 1]

    cat("   Criterion Cb (differential metabolic analysis) ...","\n")

    #sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))

    dir.create(paste0("OUTPUT-NOREVA-Criteria.Cb"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cb/Criteria.Cb-",names(normal_data[mmm]),".pdf",sep=""))
    classcolor <- NA
    for (i in 1:length(unique(obj_kmeans$cluster))){
      majority <- obj_kmeans$cluster[obj_kmeans$cluster==i]
      majority1 <- names(majority)
      majority2 <- data_kmeans[majority1,"Color"]
      majority3 <- unique(majority2)[which.max(tabulate(match(majority2, unique(majority2))))]
      majority4 <- as.character(majority3)
      classcolor[i] <- majority4
    }
    ##########################################################################################
    #################### output
    ##########################################################################################
    try(print(autoplot(obj_kmeans,
                       data = data_kmeans,
                       frame = TRUE,shape=1)+
                ggrepel::geom_text_repel(aes(label=data_kmeans$T_label),color=data_kmeans$Color,size=4,family="serif")+
                theme(legend.position="none",panel.background = element_blank(),axis.title.x=element_text(angle=0, size=12,color="black"),axis.text.x=element_text(angle=0, size=13,color="black"),axis.title.y=element_text(size=12,color="black"),axis.text.y=element_text(size=13,color="black"),panel.border = element_rect(fill='transparent', color='black'))+
                geom_point(color=data_kmeans$Color,size=3)+
                #scale_fill_manual(values = c("#BC4D70", "#00B1C1", "#55A51C")) +
                #scale_color_manual(values = c("#BC4D70", "#00B1C1", "#55A51C"))
                scale_fill_manual(values = classcolor) +
                scale_color_manual(values = classcolor)))
    dev.off()
    #####################################################################################################
    ###3. Consistency -------------------------------------------------------- #
    #####################################################################################################
    test_data <- eva_data3
    test_data[sapply(test_data, simplify = 'matrix',is.nan)] <- 0

    test_data <- test_data[order(test_data$Group),]

    number_labels <- test_data$Group
    folds <- 3
    test.fold <- list()

    test.fold1 <- sampling::strata(c("Group"),size=(as.numeric(table(test_data$Group))/3),method="srswor",data=test_data)[,2]
    test.fold[[1]] <- test.fold1

    data.2 <- test_data[-test.fold1,]

    test.fold2 <- sampling::strata(c("Group"),size=(as.numeric(table(test_data$Group))/3),method="srswor",data=data.2)[,2]
    test.fold2 <- match(row.names(data.2)[test.fold2],row.names(test_data))

    test.fold[[2]] <- test.fold2

    test.fold[[3]] <- (1:nrow(test_data))[-c(test.fold1,test.fold2)]

    tryCatch({
      DEG.list <- consistency_M(3, 80)}, error=function(e){})

    CW_value <- try(CWvalue(DEG.list,Y=(ncol(eva_data3)-1),n=length(DEG.list[[1]])))
    if(class(CW_value)=="try-error")
    { return(NA) }
    ##########################################################################################
    #################### output
    ##########################################################################################
    # Fscore[names(normal_data[mmm])]<-CW_value # output-Fscore[names(normal_data[mmm])]

    setlist3 <- DEG.list

    #sink()

    cat("   Criterion Cc (consistency in marker discovery) ...","\n")

    #sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))

    dir.create(paste0("OUTPUT-NOREVA-Criteria.Cc"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cc/Criteria.Cc-",names(normal_data[mmm]),".pdf",sep=""))
    try(print(plot(eulerr::venn(DEG.list),fills = list(fill = c("white", "white","white")),
                   labels = list(col = "black", font = 2),
                   edges = list(col = c("#800080", "#4285f4", "#fbbc05"), lwd=4),
                   quantities = TRUE)))
    dev.off()
    #####################################################################################################
    # -- 4. AUC value --------------------------------------------------------------------------------- #
    #####################################################################################################
    set.seed(3)

    X <- X_matrix[, pos_filter]
    y <- as.factor(Group)

    folds <- 5

    test.fold <- split(sample(1:length(y)), 1:folds) #ignore warning
    all.pred.tables <-  lapply(1:folds, function(i) {
      test <- test.fold[[i]]

      Xtrain <- try(X[-test, ])

      ytrain <- as.factor(y[-test])

      sm <- try(e1071::svm(Xtrain, ytrain, cost = 1000, prob = TRUE)) # some tuning may be needed

      prob.benign <- try(attr(predict(sm, X[test,], prob = TRUE), "probabilities")[, 2])

      data.frame(ytest = y[test], ypred = prob.benign) # returning this
    })

    full.pred.table <- try(do.call(rbind, all.pred.tables))
    if(class(full.pred.table)=="try-error")
    { return(NA) }

    svm_para3 <- c(1, 5, as.numeric(1))
    roc_data3 <- full.pred.table
    auc.value <- try(AUC::auc(AUC::roc(full.pred.table[, 2], full.pred.table[, 1])))
    if(class(auc.value)=="try-error")
    { return(NA) }
    ##########################################################################################
    #################### output
    ##########################################################################################
    # Fauc[names(normal_data[mmm])]<-auc.value # output-Fauc[names(normal_data[mmm])]

    ##sink()
    ##########################################################################################
    #################### output
    ##########################################################################################
    cat("   Criterion Cd (classification accuracy) ...","\n")
    #cat("\n")
    ##sink()
    dir.create(paste0("OUTPUT-NOREVA-Criteria.Cd"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cd/Criteria.Cd-",names(normal_data[mmm]),".pdf",sep=""))
    try(plot(AUC::roc(full.pred.table[, 2], full.pred.table[, 1]), col = "red"))
    dev.off() # output-pdf
    #################   save data

    ######################################################
    ###############与multi_qcs相比，新增的部分##################
    ######################################################
    # -- 5. Accuracy FC --------------------------------------------------------------------------------- #
    #################work！！！！！！！！！！！
    #################work！！！！！！！！！！！
    #################work！！！！！！！！！！！
    #################work！！！！！！！！！！！
    #################work！！！！！！！！！！！
    #################work！！！！！！！！！！！GS赋值
    path_1 <- GS
    pre_file2_1 <- readLines(path_1, n = 2)
    loc <- which.max(c(length(unlist(strsplit(pre_file2_1, ","))), length(unlist(strsplit(pre_file2_1, ";"))), length(unlist(strsplit(pre_file2_1, "\t")))))
    sep_seq <- c(",", ";", "\t")
    M_log_up <- read.csv(path_1,header=TRUE,sep=sep_seq[loc])

    if(length(intersect(colnames(M_log_up),colnames(data_q)))<2){
      stop("Criteria E cannot performed, due to without matched golden standards metabolites :\n")
    }
    #M_log_up <- read.csv("GS_time_series_QC_sim3.csv", header=TRUE)
    bias.marker <- M_log_up
    bias.marker <- M_log_up[, -1]
    rownames(bias.marker) <- M_log_up[, 1]
    colnames(bias.marker)[1] <- "Group"

    bias.marker[sapply(bias.marker, simplify = 'matrix', is.na)] <- 0.0001
    bias.marker[sapply(bias.marker, simplify = 'matrix', is.infinite)] <- 0.0001

    class_sample <- names(table(bias.marker$Group))

    bias.norm <- eva_data3

    bias.norm[sapply(bias.norm, simplify = 'matrix', is.na)] <- 0.0001
    bias.norm[sapply(bias.norm, simplify = 'matrix', is.infinite)] <- 0.0001

    class_sample2 <- names(table(bias.norm$Group))

    result<-NULL
    #i=2
    for(i in 2:(length(class_sample))){

      ##
      c1 <- bias.marker[bias.marker$Group == class_sample[i-1], -1]
      c2 <- bias.marker[bias.marker$Group == class_sample[i], -1]

      c1_mean <- apply(c1, 2, mean)
      c2_mean <- apply(c2, 2, mean)

      true_fc <- c2_mean / c1_mean
      fc_marker <- log2(true_fc) #


      names(fc_marker) <- colnames(bias.marker)[-1]

      ##
      c1_3 <- bias.norm[bias.norm$Group == class_sample2[i-1], -1]
      c2_3 <- bias.norm[bias.norm$Group == class_sample2[i], -1]


      c1_3_mean <- apply(c1_3, 2, mean)
      c2_3_mean <- apply(c2_3, 2, mean)

      if (j == 2){
        fc_norm <- c2_3_mean - c1_3_mean
      }else{
        norm_fc <- c2_3_mean / c1_3_mean
        fc_norm <- log2(norm_fc)
      }

      ##compare_fc
      mark_nor <- fc_norm[match(names(fc_marker), names(fc_norm))]
      mark_true <-  fc_marker
      bias_logfc <- mark_nor - mark_true

      if (length(unique(class_sample)) == 2){
        fc_table <- cbind(mark_true, mark_nor, bias_logfc)
        colnames(fc_table) <- c("Reference logFC", "Normalized logFC", paste(class_sample[i], "/", class_sample[i-1],sep=" "))
        #result[[i]]<-fc_table
        result<-cbind(result,fc_table[,3])

      }else if(length(unique(class_sample)) > 2){
        fc_table <- as.matrix(bias_logfc)
        colnames(fc_table) <- paste(class_sample[i], "/", class_sample[i-1],sep=" ")
        #result[[i]]<-fc_table
        result<-cbind(result,fc_table)

      }

    }

   # FCmedian <- abs(median(result,na.rm=TRUE))
   # Fcmed[names(normal_data[mmm])]<-FCmedian
    cat("   Criterion Ce (Level of Correspondence Between Normalized and Reference Data) ...","\n")
    cat("\n")

    dir.create(paste0("OUTPUT-NOREVA-Criteria.Ce"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Ce/Criteria.Ce-",names(normal_data[mmm]),".pdf",sep=""))

    #sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
    cut_col<-ggsci::pal_simpsons("springfield")(12)
    try(boxplot(result, col = cut_col[2:(length(class_sample)+1)],
                ylab = "Logarithmic Fold Change of Means",
                sub = "Differences between Two Classes"))
    abline(h = 0, col="black",lwd=0.5)

    #sink()
    dev.off()

  k.temp <- c(names(normal_data[mmm]),mean(pmad3N),accuracy,CW_value,auc.value,abs(median(result,na.rm=TRUE)))
  return(k.temp)
}

# cluster <- makeCluster(parallel::detectCores()-1, type = "SOCK") ;cluster%>% registerDoSNOW ; time = proc.time() #
# # k.test <- foreach::foreach (i = 1:length(normal_data),.packages=c("reshape")) %dopar% k.pp(i,afterqc.table,normal_data,nanmes_right,newname) # length(normal_data)
# k.test <- foreach::foreach (i = 1:50,.packages=c("reshape","ggplot2"),.errorhandling = c("pass"),.combine = "rbind") %dopar% k.pp(i,aftetable,normal_data,nanmes_right,newname) # length(normal_data)
print(proc.time()-time)
stopCluster(cluster)

save(k.test,file="./step2_data.Rdata")

k.result1 <- k.test%>%.[,-1]%>%apply(.,2,as.numeric)%>% data.table(.,id=k.test%>%.[,1]) %>% .[is.na(V1)==F]
k.result2 <- k.result1[,.(V1,V2,V3,V4,V5)] %>% as.data.frame()
colnames(k.result2) <- c("Precision","Cluster_accuracy","Reproducibility","Classification","Accuracy")
rownames(k.result2) <- k.result1$id
result2 <- k.result2

  #需要加载的包data.table,tidyverse,plyr,doSNOW,foreach,parallel

  # =========================================排名csv及排名热图的输出==============



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
  colnames(zuihou3) <- c("Overall-Rank","Criteria.Ca-Rank","Criteria.Cb-Rank","Criteria.Cc-Rank","Criteria.Cd-Rank","Criteria.Ce-Rank","Criteria.Ca-Value","Criteria.Cb-Value","Criteria.Cc-Value","Criteria.Cd-Value","Criteria.Ce-Value")
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
  colnames(data_heat) <- c("Criterion Ca: Reduction of Intragroup Variation","Criterion Cb: Differential Metabolic Analysis","Criterion Cc: Consistency in Marker Discovery","Criterion Cd: Classification Accuracy","Criterion Ce: Level of Correspondence Between Normalized and Reference Data")

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
