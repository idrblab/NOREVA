#' @title Time-course Metabolomic Study with dataset with Internal Standards (ISs).
#' @description this function enables the performance assessment of metabolomic data
#' processing for time-course dataset (with internal standards but without quality control sample)
#' using four criteria, and can scan the customized workflows selected by user and rank them based on performances.
#' @param fileName Allows the user to indicate the NAME of peak table resulted from PrepareInuputFiles() (default = null).
#' @param selectedMethods Allows the user to indicate the NAME of the file containing the
#' customized workflows selected by user. The file should be in a .csv format,
#' and the exemplar files are provided in the NOREVA R package and available for
#' download at https://idrblab.org/noreva/NOREVA_exampledata.zip).
#' @param IS Allows the user to indicate the column number(s) where the internal standard(s) locate (default = null).
#' If there is only one internal standard (IS), the column number of this IS should be listed.
#' If there are multiple ISs, the column numbers of all ISs should be listed and separated using comma.
#' For example, the value of argument IS that is set to “2,6,8,n” indicates that the metabolites in the 3rd, 7th, 9th, and (n+1)th columns of your input peak table should be considered to be the IS metabolites.
#' @import DiffCorr affy vsn DT
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
#' @usage nortimecourseispart(fileName, IS, selectedMethods)
#' @export nortimecourseispart
#' @examples
#' library(NOREVA)
#' \donttest{timec_is_data <- PrepareInuputFiles(dataformat = 1,
#' rawdata = "Timecourse_with_IS.csv")}
#' \donttest{nortimecourseispart(fileName = timec_is_data,
#' , IS = "4,5", selectedMethods = "selectedmethods.csv")}


nortimecourseispart <- function(fileName, IS, selectedMethods){

  cat("NOREVA is Running ...","\n")
  cat("\n")
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


  norm_no <- "NON"
  #-----------------------------------------------------------------

  ###################################################Step-2 调用数据
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

  path_1 <- selectedMethods
  pre_file2_1 <- readLines(path_1, n = 2)
  loc <- which.max(c(length(unlist(strsplit(pre_file2_1, ","))), length(unlist(strsplit(pre_file2_1, ";"))), length(unlist(strsplit(pre_file2_1, "\t")))))
  sep_seq <- c(",", ";", "\t")
  data_p <- read.csv(path_1,header=TRUE,sep=sep_seq[loc])
  #data_p <- selectFile

  for(p in 1:nrow(data_p)){
    #p = 2
    datap1 <- data_p[p,2]
    datap2 <- as.character(datap1)
    datap3 <- unlist(strsplit(datap2,"+",fixed= T))

    impt <- datap3[1]; trsf <- datap3[2]; nmals <- datap3[3]

    message(sprintf("impt:%s; ", impt),sprintf("trsf:%s; ", trsf),sprintf("nmals:%s; ", nmals))


  #for (i in as.numeric(impt)){
    for (i in 1:4){

    after.table <- NULL
    imput_m <- try(imput(filter_data2,i))
    if(class(imput_m)=="try-error")
    { next }
    train_data<- t(imput_m)

    #for (j in as.numeric(trsf)){
      for (j in 1:3){

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

        normal_data[[paste(im_nam[i], norm_no, t_nam[j], paste("[",paste(n_norm_IS[k]),"]", sep=""), sep="+")]] <- eva_data3
        eva_data3 <- NULL
        aftetable[[paste(im_nam[i],t_nam[j],n_norm_IS[k],sep="+")]] <- after.table

      }
    }
  }
}
  #length(normal_data)
  nanmes_right<-names(normal_data)
  save(normal_data,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
  load("./OUTPUT-NOREVA-All.Normalized.Data.Rdata")

  Fpmad<-list()
  Fpurity<-list()
  Fscore<-list()
  Fauc<-list()

  for (mmm in 1:length(normal_data)){

    name <- nanmes_right

    ###Fpmad-------------------------------------------------------------------------
    n_data <- as.data.frame(normal_data[mmm],col.names=NULL)
    n_data <- as.matrix(n_data)
    if(sum(is.na(n_data))<length(n_data)/3){
      eva_data3<-as.data.frame(normal_data[mmm],col.names=NULL)
    }else{next}

    eva_data3_f<-eva_data3

    eva_data3 <- eva_data3_f

    ###

    pmad3N.log <- eva_data3
    pmad3N <- try(PMAD(pmad3N.log))
    if(class(pmad3N)=="try-error")
    { next }

    Fpmad[names(normal_data[mmm])]<-mean(pmad3N)

    names(after.table)[1] <- "Group"

    pmad3R.log <- after.table[,-1]
    pmad3R <- try(PMAD(pmad3R.log))
    if(class(pmad3R)=="try-error")
    { next }

    C3 <- cbind(pmad3R, pmad3N); colnames(C3) <- c("Before", "After")

    cat(paste("Assessing Method" , paste(mmm,":",sep=""), name[mmm]),"\n")

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
    # 2. Fpurity-------------------------------------------------------------------------

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
    { next }

    marker_ranking <- cbind(gnames, out1$HotellingT2)

    DEG <- marker_ranking[order(as.numeric(marker_ranking[,2]), decreasing=T),]

    data_kmeans <- as.data.frame(eva_data3_f[, c("Group",DEG[(1:50),1])])

    clusters <- length(unique(data_kmeans[, 1]))
    obj_kmeans <- try(kmeans(data_kmeans[,-1], centers = clusters, nstart = 20,iter.max = 100))
    if(class(obj_kmeans)=="try-error")
    { next }

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

    accuracy<-purity(result,label)
    Fpurity[names(normal_data[mmm])]<-accuracy

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

    # 3. Consistency -------------------------------------------------------- #

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

    DEG.list<-try(consistency(3, 20))
    if(class(DEG.list)=="try-error")
    { next }

    CW_value <- try(CWvalue(DEG.list,Y=(ncol(eva_data3)-1),n=length(DEG.list[[1]])))
    if(class(CW_value)=="try-error")
    { next }

    Fscore[names(normal_data[mmm])]<-CW_value

    sink()

    cat("   Criterion Cc (consistency in marker discovery) ...","\n")


    sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
    dir.create(paste0("OUTPUT-NOREVA-Criteria.Cc"))
    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cc/Criteria.Cc-",names(normal_data[mmm]),".pdf",sep=""))
    try(print(plot(eulerr::venn(DEG.list),fills = list(fill = c("white", "white","white")),
                   labels = list(col = "black", font = 2),
                   edges = list(col = c("#800080", "#4285f4", "#fbbc05"), lwd=4),
                   quantities = TRUE)))
    dev.off()
    # -- 4. AUC value --------------------------------------------------------------------------------- #

    DEG <- marker_ranking[order(as.numeric(marker_ranking[,2]), decreasing=T),]

    set.seed(3)

    # NB ROC PLOT will change for each new random noise component (jitter)
    X_matrix <- eva_data3[,-1]
    y_label <- as.factor(eva_data3[,1])
    X_matrix[sapply(X_matrix, simplify = 'matrix', is.na)] <- 0

    data_multiROC <- cbind(y_label, X_matrix)
    data_multiROC <- cbind(data_multiROC[,-1], paste("Label_", data_multiROC[, 1], sep = ""))
    colnames(data_multiROC)[ncol(data_multiROC)] <- "Label"

    X_matrix <- data_multiROC
    y_label <- data_multiROC[,"Label"]

    x <- as.data.frame(X_matrix[, c(DEG[(1:20),1], "Label")])
    y <- y_label

    #cross validation
    y<- as.factor(x[, length(x)])
    folds <- 5
    test.fold <- split(sample(1:length(y)), 1:folds) #ignore warning

    for (mm in 1:5) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]

      svmmodel <- try(svm(Label ~ ., data = train_df, probability = TRUE))
      if(class(svmmodel)=="try-error")
      { next }

      svm_pred <- try(predict(svmmodel, test_df,probability = TRUE))
      if(class(svm_pred)=="try-error")
      { next }

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
    roc_res <- try(multi_roc(final_df_m_t, force_diag = T))
    if(class(roc_res)=="try-error")
    { next }

    plot_roc_df <- try(plot_roc_data(roc_res))
    if(class(plot_roc_df)=="try-error")
    { next }

    plot_roc_df_mic<-subset(plot_roc_df, plot_roc_df$Group=="Micro")

    AUC_mic<-round(unique(plot_roc_df_mic$AUC),3)

    Fauc[names(normal_data[mmm])]<-AUC_mic

    sink()

    cat("   Criterion Cd (classification accuracy) ...","\n")
    cat("\n")
    dir.create(paste0("OUTPUT-NOREVA-Criteria.Cd"))

    pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cd/Criteria.Cd-",names(normal_data[mmm]),".pdf",sep=""))
    plot(x=c(0,1),y=c(0,1),col="lightgrey",pch=16,bg="yellow",type = 'l',xlim=c(0,1),ylim=c(0,1),lwd=1,xlab="1-Specificity",ylab="Sensitivity")
    lines(x=1-plot_roc_df_mic$Specificity,plot_roc_df_mic$Sensitivity,col="red",pch=16,bg="yellow",xlim=c(0,1),ylim=c(0,1),lwd=2,xlab="WEEK",ylab="STUDE")

    dev.off()
  }

  result<-dplyr::bind_rows("Precision"=unlist(Fpmad),"Cluster_accuracy"=unlist(Fpurity),"Reproducibility"=unlist(Fscore),"Classification"=unlist(Fauc),.id = "id")
  result1<-t(result)
  colnames(result1)<-result1["id",]
  result2<-result1[-1,]
  result2<-data.frame(result2,check.names=FALSE)

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
  cat("\n")
  cat("*************************************************************************","\n")
  cat("Congratulations! Assessment Successfully Completed!","\n")
  cat("Thanks for Using NOREVA. Wish to See You Soon ...","\n")
  cat("*************************************************************************","\n")
  cat("\n")
  #return(rank_result)
}
