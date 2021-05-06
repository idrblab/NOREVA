#' @title Process the dataset of time-course metabolomic study.
#' @description based on a particular processing workflow
#' (especially the one identified as well-performing for the studied metabolomic time-course dataset),
#' this function outputs the resulting levels of all metabolites among all samples after the
#' data processing based on that workflow.
#' @param datatype Input the number of data type.
#' If set 1, the dataset of time-course metabolomic study without QCSs and ISs.
#' If set 2, the dataset of time-course metabolomic study with quality control samples (QCSs).
#' If set 3, the dataset of time-course metabolomic study with dataset with internal standards (ISs).
#' @param fileName Please find the detail information of the file format from those six sample
#' files in the working directory (in github) “idrblab/NOREVA/data”
#' @param IS Input the Column of Internal Standards. For example,
#' the replacement of IS to 2,6,9,n indicates that the metabolites in the 2st, 6th, 9th, and nth columns of in
#' your input dataset Input-Dataset.csv should be considered as the ISs or quality control metabolites.
#' If there is only one IS, the column number of this IS should be listed
#' If there are multiple ISs, the column number of all ISs should be listed and
#' separated by comma (,)
#' @param impt Input the name of imputation methods.
#' If set 1, method of column mean imputation.
#' If set 2, method of column median imputation.
#' If set 3, method of half of the minimum positive value.
#' If set 4, method of KNN imputation.
#' @param qcsn Input the name of qc sample correction methods.
#' If set 1, method of NWE (Nadaraya-Watson estimator).
#' If set 2, method of LLR (local linear regression).
#' If set 3, method of LPF (local polynomial fits).
#' @param trsf Input the name of transformation methods.
#' If set 1, method of cube root transformation.
#' If set 2, method of log transformation.
#' If set 3, none transformation method.
#' @param nmal Allows the users to specify the NAME of the normalization method (default = null)
#' @param nmal2 Allows the users to specify the NAME of the normalization method (default = null)
#' @param nmals Allows the users to specify the NAME of the IS-based normalization method (default = null)
#' @import DiffCorr affy vsn DT
#' @import e1071 AUC impute MetNorm
#' @import ggsci timecourse multiROC dummies
#' @import ggplot2 ggord ggfortify usethis
#' @import ggrepel ggpubr sampling crmn
#' @rawNamespace import(limma, except=.__C__LargeDataObject)
#' @rawNamespace import(ropls, except=plot)
#' @importFrom grDevices dev.off png rainbow rgb colorRampPalette pdf
#' @importFrom graphics abline close.screen legend mtext par points screen split.screen symbols text title
#' @importFrom stats anova as.formula cor dnorm kmeans lm loess loess.control mad median model.matrix na.omit pf pnorm qnorm qt quantile rnorm runif sd var
#' @importFrom utils combn read.csv write.csv write.table
#' @examples nordata <- nortimecoursematrix(datatype = 1,
#' fileName = timec_non_data, impt = 1, trsf = 1, nmal = 1, nmal2 = 1)
#' @export nortimecoursematrix

nortimecoursematrix <- function(datatype, fileName, IS, impt=NULL, qcsn=NULL, trsf=NULL, nmal=NULL, nmal2=NULL, nmals=NULL){

  cat("NOREVA is Running ...","\n")
  cat("\n")

  #imputation---------------------------------------------------------------------------------
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
    "Mean",
    "Median",
    "Half Minimum",
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
    "Cube_root",
    "Log",
    "None"
  )

  q_nam<-c(
    "NWE",
    "LLR",
    "LPF"
  )

  #normalization---------------------------------------------------------------------------------
  norm_sam<-function(train_data_t,n){
    matrix<-switch(
      n,
      train_data_t,#1
      PQN(train_data_t),#2
      fastlo(train_data_t),#3
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

  #normalization for with IS---------------------------------------------------------------------------------
  norm_IS<-function(train_data_nor,n){
    matrix<-switch(
      n,
      t(SIS(train_data_nor, as.numeric(nomis_name))),#1
      t(NOMIS(train_data_nor, as.numeric(nomis_name))),#2
      t(CCMN(train_data_nor, as.numeric(nomis_name))),#3
      t(RUVRand(train_data_nor, as.numeric(nomis_name)))#4
    )
    return(matrix)
  }
  n_norm_IS<-c(
    "SIS",
    "NOMIS",
    "CCMN",
    "RUV-random"
  )

  #normalization for with IS---------------------------------------------------------------------------------
  #-----------------------------------------------------------------
  ###################################################Step-2 调用数据
  if (as.numeric(datatype) == 1){
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

    impt <- impt; trsf <- trsf; nmal <- nmal; nmal2 <- nmal2

    #message(sprintf("impt:%s; ", impt),sprintf("trsf:%s; ", trsf),sprintf("nmal:%s; ", nmal),sprintf("nmal2:%s\t", nmal2))


    for (i in as.numeric(impt)){
      #for (i in 1:7){

      afterimpute.table <- NULL
      imput_m <- try(imput(filter_data2,i))
      if(class(imput_m)=="try-error")
      { next }

      imputed_data <- cbind(data4, imput_m)
      afterimpute.table <- imputed_data
      getLabel <-as.character(data_q[, 2])
      data1 <- afterimpute.table
      train_data_t <- t(data1[, -(1:2)])

      for (j in as.numeric(trsf)){
        #for (j in 1:3){

        train_data_Transformation3<-try(trans(train_data_t,j))
        if(class(train_data_Transformation3)=="try-error")
        { next }
        #train_data_t <- train_data_Transformation3
        train_data_Transformation3[is.infinite(data.matrix(train_data_Transformation3))]<-NA
        sampleLabel <- as.character(afterimpute.table[, 2])

        imputed_data <- cbind(data4, t(train_data_Transformation3))
        after.table <- imputed_data

        ####
        if (nmal == 1|nmal == 14){
          for(k in as.numeric(nmal)){

            train_data_metabolite3<-try(norm_sam(train_data_Transformation3,k))
            if(class(train_data_metabolite3)=="try-error")
            { next }

            normalized_data3 <- try(t(train_data_metabolite3))
            if(class(normalized_data3)=="try-error")
            { next }

            eva.data3 <- cbind(afterimpute.table[, 1:2], normalized_data3)
            eva.data3 <- eva.data3[, -1]
            eva.data3 <- as.data.frame(eva.data3)
            rownames(eva.data3) <- afterimpute.table[, 1]
            colnames(eva.data3)[1] <- "Group"
            eva_data3<-eva.data3

            normal_data[[paste(im_nam[i],t_nam[j],n_nam[k],sep="+")]] <- eva_data3

            #eva_data3 <- NULL
            aftetable[[paste(im_nam[i],t_nam[j],n_nam[k],sep="+")]] <- after.table

          }
        }else {
          for ( k in as.numeric(nmal)){
            train_data_Preprocess3 <-try(norm_sam(train_data_Transformation3,k))
            if(class(train_data_Preprocess3)=="try-error")
            { next }

            for(h in as.numeric(nmal2)){

              train_data_metabolite3<-try(norm_sam(train_data_Preprocess3,h))
              if(class(train_data_metabolite3)=="try-error")
              { next }

              normalized_data3 <- try(t(train_data_metabolite3))
              if(class(normalized_data3)=="try-error")
              { next }

              eva.data3 <- cbind(afterimpute.table[, 1:2], normalized_data3)
              eva.data3 <- eva.data3[, -1]
              eva.data3 <- as.data.frame(eva.data3)
              rownames(eva.data3) <- afterimpute.table[, 1]
              colnames(eva.data3)[1] <- "Group"
              eva_data3<-eva.data3

              normal_data[[paste(im_nam[i],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- eva_data3

              #eva_data3 <- NULL

              aftetable[[paste(im_nam[i],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- after.table

            }
          }
        }
      }
    }
  }else if(as.numeric(datatype) == 2){

    data_q <- fileName

    sampleData <- as.data.frame(data_q[, -(1:5)])
    sampleLabel_00 <- as.character(data_q[, 1])
    sampleList <- as.data.frame(data_q[, 1:4])
    names(sampleList) <- c("sample", "batch", "class", "order")
    time <- as.character(data_q[-which(is.na(data_q[,3])), 5])

    data2 <- sampleData
    data2_QC <- data_q[which(is.na(data_q[,3])),-(1:5)]
    data2_QC[data2_QC == 0] <- NA

    #filtering-----------------------------------------------------------------------------
    col_f <- apply(data2_QC, 2, function(x) length(which(is.na(x)))/length(x))

    if (length(which(col_f > 0.2))==0){
      data2_f <- data2
    }else {
      data2_f <- data2[, -which(col_f > 0.2)]
    }


    col_r <- apply(data2_QC, 2, function(x) sd(x)/mean(x))

    if (length(which(col_r > 0.3))==0){
      filter_data2 <- data2_f
    }else {
      filter_data2 <- data2_f[, -which(col_r > 0.3)]
    }


    #---------------------------------------------------------------------------------

    normal_data <- list()
    aftetable<-list()
    eva_data3 <- NULL

    impt <- impt; trsf <- trsf; nmal <- nmal; nmal2 <- nmal2

    #message(sprintf("impt:%s; ", impt),sprintf("trsf:%s; ", trsf),sprintf("nmal:%s; ", nmal),sprintf("nmal2:%s\t", nmal2))

    for (i in as.numeric(impt)){
      #for (i in 1:7){

      afterqc.table<-NULL
      imput_m2 <- try(imput(filter_data2,i))
      if(class(imput_m2)=="try-error")
      { next }

      # running the QC processing.
      imput_m2_t<-t(imput_m2)
      sampleData_rev <- as.data.frame(cbind(rownames(imput_m2_t), imput_m2_t))
      colnames(sampleData_rev)<- c("name",sampleLabel_00)
      samPeno <- sampleList
      samFile <- sampleData_rev

      for (q in as.numeric(qcsn)){

      degree <- q-1

      afterqc_table <- try(shiftCor(samPeno, samFile, Frule = 0.8, QCspan = 0.75, degree = degree))

      if(class(afterqc_table)=="try-error")
      { next }

      afterqc_table <- cbind(afterqc_table[,1:2],time, afterqc_table[,-(1:2)])
      afterqc_table<-afterqc_table[,-2]
      afterqc.table<-afterqc_table[order(afterqc_table[,2]),]
      getLabel <-as.character(afterqc.table[, 2])
      data1 <- afterqc.table
      train_data_t <- t(data1[, -(1:2)])

      for (j in as.numeric(trsf)){
        #for (j in 1:3){

        train_data_Transformation<-try(trans(train_data_t,j))
        if(class(train_data_Transformation)=="try-error")
        { next }
        sampleLabel <- as.character(afterqc.table[, 2])

        imputed_data <- cbind(afterqc_table[,1:2],t(train_data_Transformation))

        afterqc.trans <- imputed_data

        ####
        if (nmal == 1|nmal == 14){

          for(k in as.numeric(nmal)){
            train_data_metabolite3<-try(norm_sam(train_data_Transformation,k))
            if(class(train_data_metabolite3)=="try-error")
            { next }

            normalized_data3 <- try(t(train_data_metabolite3))
            if(class(normalized_data3)=="try-error")
            { next }

            eva.data3 <- cbind(afterqc.table[, 1:2], normalized_data3)
            eva.data3 <- eva.data3[,-1]
            eva.data3 <- as.data.frame(eva.data3)
            rownames(eva.data3) <- afterqc.table[, 1]
            colnames(eva.data3)[1] <- "Group"
            eva_data3<-eva.data3

            normal_data[[paste(im_nam[i],q_nam[q],t_nam[j],n_nam[k],sep="+")]] <- eva_data3
            #eva_data3 <- NULL
            aftetable[[paste(im_nam[i],q_nam[q],t_nam[j],n_nam[k],sep="+")]] <- afterqc.trans}}else {

              for ( k in as.numeric(nmal)){
                train_data_Preprocess3 <-try(norm_sam(train_data_Transformation,k))
                if(class(train_data_Preprocess3)=="try-error")
                { next }
                for(h in as.numeric(nmal2)){
                  train_data_metabolite3<-try(norm_sam(train_data_Preprocess3,h))
                  if(class(train_data_metabolite3)=="try-error")
                  { next }
                  #if (train_data_metabolite3==NULL){
                  #  normal_data[[paste(im_nam[i], t_nam[j], n_nam[k],n_nam[h],sep="+")]] <- "abc"
                  #}else{
                  normalized_data3 <- try(t(train_data_metabolite3))
                  if(class(normalized_data3)=="try-error")
                  { next }

                  eva.data3 <- cbind(afterqc.table[, 1:2], normalized_data3)
                  eva.data3 <- eva.data3[,-1]
                  eva.data3 <- as.data.frame(eva.data3)
                  rownames(eva.data3) <- afterqc.table[, 1]
                  colnames(eva.data3)[1] <- "Group"
                  eva_data3<-eva.data3

                  normal_data[[paste(im_nam[i],q_nam[q],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- eva_data3
                  #eva_data3 <- NULL
                  aftetable[[paste(im_nam[i],q_nam[q],t_nam[j],n_nam[k],n_nam[h],sep="+")]] <- afterqc.trans
                }
              }
            }
      }
      }
    }
  }else if(as.numeric(datatype) == 3){
    message("Please check IS should be a series of natural numbers separated by comma.")
    internal_standard <- IS
    data_q <- fileName

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

    aftetable<-list()
    normal_data <- list()
    train_data_Preprocess3 <- NULL

    impt <- impt; trsf <- trsf; nmals <- nmals

    #message(sprintf("impt:%s; ", impt),sprintf("trsf:%s; ", trsf),sprintf("nmals:%s; ", nmals))

    for (i in as.numeric(impt)){
      #for (i in 1:7){

      after.table <- NULL
      imput_m <- try(imput(filter_data2,i))
      if(class(imput_m)=="try-error")
      { next }
      train_data<- t(imput_m)

      for (j in as.numeric(trsf)){
        #for (j in 1:3){

        train_data_Transformation3<-try(trans(train_data,j))
        if(class(train_data_Transformation3)=="try-error")
        { next }

        train_data_Transformation3[is.infinite(data.matrix(train_data_Transformation3))]<-NA
        after.table <- cbind(data4, t(train_data_Transformation3))

        train_data_tr<-after.table[,-1]
        sampleLabel <- as.character(after.table[, 2])

        for ( k in as.numeric(nmals)){
          #Data format:(1)Matrix: row is samples, column is metabolites.(The first column is the binary labels.)
          ###            (2)nc is the column order of QC metabolites or IS.

          is_name <- as.numeric(unlist(strsplit(as.character(internal_standard),",")))

          nomis_name <- is_name
          #nomis_name=c(2,3,4)###################################

          train_data_Preprocess3 <-try(norm_IS(train_data_tr,k))
          if(class(train_data_Preprocess3)=="try-error")
          { next }
          normalized_data3 <- train_data_Preprocess3

          eva.data3 <- cbind(after.table[, 1:2], t(normalized_data3))
          eva.data3 <- eva.data3[, -1]
          eva.data3 <- as.data.frame(eva.data3)
          rownames(eva.data3) <- after.table[, 1]
          colnames(eva.data3)[1] <- "Group"
          eva_data3<-eva.data3

          normal_data[[paste(im_nam[i],t_nam[j],n_norm_IS[k],sep="+")]] <- eva_data3
          #eva_data3 <- NULL
          aftetable[[paste(im_nam[i],t_nam[j],n_norm_IS[k],sep="+")]] <- after.table

        }
      }
    }
  }
  #nanmes_right<-names(normal_data)
  #write.csv(eva_data3,file = "./OUTPUT-NOREVA-Normalized.Data.csv")
  cat("\n")
  cat("*************************************************************************","\n")
  cat("Congratulations! Assessment Successfully Completed!","\n")
  cat("Thanks for Using NOREVA. Wish to See You Soon ...","\n")
  cat("*************************************************************************","\n")
  cat("\n")
  return(eva_data3)
}
