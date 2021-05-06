#' @title Multi-class (N>1) Metabolomic Study with dataset with Internal Standards (ISs) and the corresponding data of golden standards for performance evaluation using Criterion e.
#' @description this function enables the performance assessment of metabolomic
#' data processing for multi-class dataset (with internal standards but without
#' quality control sample) using five criteria, and can scan thousands of processing
#' workflows and rank them based on their performances.
#' @param fileName Allows the user to indicate the NAME of peak table resulted from PrepareInuputFiles() (default = null).
#' @param IS Allows the user to indicate the column number(s) where the internal standard(s) locate (default = null)
#' If there is only one internal standard (IS), the column number of this IS should be listed
#' If there are multiple ISs, the column numbers of all ISs should be listed and separated using comma
#' For example, the value of argument IS that is set to “2,6,8,n” indicates that the metabolites in the 3rd, 7th, 9th, and (n+1)th columns of your input peak table should be considered to be the IS metabolites.
#' @param GS Allows the user to indicate the name of the file that contains the spike-in compounds (default = null).
#' The file should be in a .csv format, which provides the concentrations of spike-in compounds.
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
#' @usage normulticlassisallgs(fileName, IS, GS)
#' @export normulticlassisallgs
#' @examples
#' library(NOREVA)
#' \donttest{multi_is_data <- PrepareInuputFiles(dataformat = 1,
#' rawdata = "Multiclass_with_IS.csv")}
#' \donttest{normulticlassisallgs(fileName = multi_is_data,
#' GS = "Multiclass_with_IS_GoldenStandard.csv", IS = "3,4,5")}

normulticlassisallgs <- function(fileName, IS, GS){

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
  cat("STEP 2: The criteira selected by users for this assessment","\n")
  cat("Criterion Ca: Reduction of Intragroup Variation", "\n")
  cat("Criterion Cb: Differential Metabolic Analysis", "\n")
  cat("Criterion Cc: Consistency in Marker Discovery", "\n")
  cat("Criterion Cd: Classification Accuracy", "\n")
  cat("Criterion Ce: Level of correspondence between normalized and reference data", "\n")
  cat("\n")

  cat("NOREVA is Running ...","\n")
  cat("\n")

  #To ignore the warnings during usage
  options(show.error.messages=FALSE,echo=FALSE,keep.source.pkg=TRUE)
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
  "NOM",
  "CCM",
  "RUV"
)

norm_no <- "NON"
#-----------------------------------------------------------------

###################################################Step-2 调用数据
internal_standard <- IS
data_q <- fileName

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

    for (k in 1:4){
      #Data format:(1)Matrix: row is samples, column is metabolites.(The first column is the binary labels.)
      ###            (2)nc is the column order of QC metabolites or IS.

      #nomis_name=c(2,3,4)###################################

      is_name <- as.numeric(unlist(strsplit(as.character(internal_standard),",")))

      nomis_name <- is_name
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

#length(normal_data)
nanmes_right<-names(normal_data)
save(normal_data,file="./OUTPUT-NOREVA-All.Normalized.Data.Rdata")
load("./OUTPUT-NOREVA-All.Normalized.Data.Rdata")

Fpmad<-list()
Fpurity<-list()
Fscore<-list()
Fauc<-list()
Fcmed <- list()

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
  sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
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
  pos_filter <- try(OPLSDA_test(X_matrix, Group, cutoff = 0.8))
  if(class(pos_filter)=="try-error")
  { next }

  DEG <- cbind(Group,X_matrix[, pos_filter])
  data_kmeans<-DEG

  clusters <- length(unique(data_kmeans[, 1]))
  obj_kmeans <- try(kmeans(data_kmeans[,-1], centers = clusters, nstart = 10,iter.max = 6))
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

  accuracy<-try(purity(result,label))
  if(class(accuracy)=="try-error")
  { next }

  Fpurity[names(normal_data[mmm])]<-accuracy

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
  test_data[sapply(test_data, simplify = 'matrix',is.nan)] <- 0

  test_data <- test_data[order(test_data$Group),]

  number_labels <- test_data$Group
  folds <- 3
  test.fold <- list()

  test.fold1 <- strata(c("Group"),size=(as.numeric(table(test_data$Group))/3),method="srswor",data=test_data)[,2]
  test.fold[[1]] <- test.fold1

  data.2 <- test_data[-test.fold1,]

  test.fold2 <- strata(c("Group"),size=(as.numeric(table(test_data$Group))/3),method="srswor",data=data.2)[,2]
  test.fold2 <- match(row.names(data.2)[test.fold2],row.names(test_data))

  test.fold[[2]] <- test.fold2

  test.fold[[3]] <- (1:nrow(test_data))[-c(test.fold1,test.fold2)]

  DEG.list <- try(consistency_M(3, 70))
  if(class(DEG.list)=="try-error")
  { next }

  CW_value <- try(CWvalue(DEG.list,Y=(ncol(eva_data3)-1),n=length(DEG.list[[1]])))
  if(class(CW_value)=="try-error")
  { next }

  Fscore[names(normal_data[mmm])]<-CW_value

  setlist3 <- DEG.list
  #OLlist3 <- try(overLapper(setlist = setlist3, sep="_", type = "vennsets"))
  #if(class(OLlist3)=="try-error")
  #{ next }
  sink()
  #counts <- list(sapply(OLlist3$Venn_List, length), sapply(OLlist3$Venn_List, length))
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

  set.seed(3)

  X <- X_matrix[, pos_filter]
  y <- as.factor(Group)

  # cross-validated SVM-probability plot
  folds <- 5

  test.fold <- split(sample(1:length(y)), 1:folds) #ignore warning
  all.pred.tables <-  lapply(1:folds, function(i) {
    test <- test.fold[[i]]
    Xtrain <- X[-test, ]
    ytrain <- as.factor(y[-test])

    sm <- try(svm(Xtrain, ytrain, cost = as.numeric(100),probability = TRUE)) # some tuning may be needed

    prob.benign <- attr(predict(sm, X[test,], prob = TRUE), "probabilities")[, 2]
    data.frame(ytest = y[test], ypred = prob.benign) # returning this
  })

  full.pred.table <- try(do.call(rbind, all.pred.tables))
  if(class(full.pred.table)=="try-error")
  { next }

  svm_para3 <- c(1, 5, as.numeric(1))

  roc_data3 <- full.pred.table

  auc.value <- try(auc(roc(full.pred.table[, 2], full.pred.table[, 1])))
  if(class(auc.value)=="try-error")
  { next }

  Fauc[names(normal_data[mmm])]<-auc.value
  sink()
  cat("   Criterion Cd (classification accuracy) ...","\n")
  #cat("\n")
  dir.create(paste0("OUTPUT-NOREVA-Criteria.Cd"))
  pdf(file=paste("./OUTPUT-NOREVA-Criteria.Cd/Criteria.Cd-",names(normal_data[mmm]),".pdf",sep=""))
  try(plot(roc(full.pred.table[, 2], full.pred.table[, 1]), col = "red"))
  dev.off()

  # -- 5. Accuracy FC --------------------------------------------------------------------------------- #
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

  FCmedian <- abs(median(result,na.rm=TRUE))
  Fcmed[names(normal_data[mmm])]<-FCmedian
  cat("   Criterion Ce (Level of Correspondence Between Normalized and Reference Data) ...","\n")
  cat("\n")

  dir.create(paste0("OUTPUT-NOREVA-Criteria.Ce"))
  pdf(file=paste("./OUTPUT-NOREVA-Criteria.Ce/Criteria.Ce-",names(normal_data[mmm]),".pdf",sep=""))

  sink(file=paste("OUTPUT-NOREVA-Record",".txt",sep=""))
  cut_col<-pal_simpsons("springfield")(12)
  try(boxplot(result, col = cut_col[2:(length(class_sample)+1)],
              ylab = "Logarithmic Fold Change of Means",
              sub = "Differences between Two Classes"))
  abline(h = 0, col="black",lwd=0.5)

  sink()
  dev.off()
  #cat("   Criterion Ce (spiked accuracy) ...","\n")
  #cat("\n")
}

result<-dplyr::bind_rows("Precision"=unlist(Fpmad),"Cluster_accuracy"=unlist(Fpurity),"Reproducibility"=unlist(Fscore),"Classification"=unlist(Fauc),"Accuracy"=unlist(Fcmed),.id = "id")
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
cat("\n")
cat("*************************************************************************","\n")
cat("Congratulations! Assessment Successfully Completed!","\n")
cat("Thanks for Using NOREVA. Wish to See You Soon ...","\n")
cat("*************************************************************************","\n")
cat("\n")
#return(rank_result)
}

