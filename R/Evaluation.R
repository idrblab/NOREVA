
###==========================================================================###
### Other evaluation indexes, such as PEV, PMAD, PCV, and so on.
###==========================================================================###

#calculate PEV
PEV <- function(data0){
    data0<-as.matrix(data0)
    data0<-data0[order(data0[,1]),]
    label <- as.factor(data0[,1])
    data <- data0[,-1]
    data <- t(as.matrix(data))
    x<-levels(label)[1]
    z<-1
    y<-1
    flag<-1
    count<-0
    varmem<-vector()
    tempvar<-vector()
    nonmissingmat<-vector()
    for(i in 1:length(label))
    {
        if(x!=label[i] || i==length(label))
        {
            y<-i-1
            if(i==length(label))
            {
                y<-i
            }
            if(flag==1)
            {
                count<-count+1
                nonmissingmat<-(apply(data[,z:y],1,function(x) {((sum(!is.na(x))))}))-1
                tempvar<-nonmissingmat*apply(data[,z:y],1,function(x) {var(x,na.rm=TRUE)})
            }
            if(flag==2)
            {
                count<-count+1
                nonmissingmat<-(apply(data[,z:y],1,function(x) {((sum(!is.na(x))))}))-1
                tempvar<-nonmissingmat*apply(data[,z:y],1,function(x) {var(x,na.rm=TRUE)})
            }
            varmem<-c(varmem,((sum(tempvar,na.rm=T))/(sum(nonmissingmat,na.rm=T))))
            z<-i
            x<-label[i]
            flag=2;
        }
    }
    avgvarmem<-varmem
    names(avgvarmem)<-levels(label)
    return(avgvarmem)
}
#calculate PMAD
PMAD <- function(data0){
    data0<-as.matrix(data0)
    data0<-data0[order(data0[,1]),]
    label <- as.factor(data0[,1])
    data <- data0[,-1]
    data <- t(as.matrix(data))
    data <- apply(data, 1:2, as.numeric)
    x<-levels(label)[1]
    z<-1
    y<-1
    flag<-1
    count<-0
    madmem<-matrix(nrow=nrow(data),ncol=length(levels(as.factor(unlist(label)))),byrow=T)
    for(i in 1:length(label))
    {
        if(x!=label[i] || i==length(label))
        {
            y<-i-1
            if(i==length(label))
            {
                y<-i
            }
            if(flag==1)
            {
                count<-count+1
                madmem[,count]<-apply(data[,z:y],1,function(x) {mad(x,na.rm=T)})
            }
            if(flag==2)
            {
                count<-count+1
                madmem[,count]<-apply(data[,z:y],1,function(x) {mad(x,na.rm=T)})
            }
            z<-i
            x<-label[i]
            flag=2;
        }
    }
    avgmadmem<-apply(madmem,2,mean,na.rm=T)
    names(avgmadmem)<-levels(label)
    return(avgmadmem)
}

#' @importFrom RcmdrMisc numSummary
#calculate PCV
PCV <- function(data0){
    data0<-as.matrix(data0)
    data0<-data0[order(data0[,1]),]
    label <- as.factor(data0[,1])
    data <- data0[,-1]
    data <- t(as.matrix(data))
    data <- apply(data, 1:2, as.numeric)

    tempcvmat<-matrix(nrow=nrow(data),ncol=length(levels(as.factor(unlist(label)))),byrow=T)
    for(i in 1:nrow(data))
    {
        tempcv<-numSummary(data[i,],statistics=c("cv"),groups=unlist(label))
        tempcvmat[i,]<-tempcv$table
    }
    temcvmatsum<-apply(tempcvmat,2,mean,na.rm=T)
    avgcvmem <-((temcvmatsum*100))
    names(avgcvmem)<-levels(label)
    return(avgcvmem)
}



  purity<-function(result,label){
    total_num<-length(label)
    cluster_counter = unique(result)
    original_counter=unique(label)
    t<-NULL

    for (k in cluster_counter){
      p_k<-NULL
      for (j in original_counter){
        count<-0
        for (i in 1:length(result)){
          if (result[i]==k && label[i]==j){
            count<-count+1
          }

        }
        p_k<-c(p_k,count)
      }
      temp_t<-max(p_k)
      t<-c(t,temp_t)
    }
    res<-sum(t)/total_num
    return(res)

  }



  CWvalue <- function(all_genelist=all_genelist,Y,n){
    num <- 3
    partial_genelist <- vector("list", length(all_genelist))
    for (j in 1:length(all_genelist)) {
      partial_genelist[[j]] <- all_genelist[[j]][1:n]
    }

    partial_genelist <- as.list(partial_genelist)
    partial_genelist <- table(unlist(partial_genelist))

    Ff_sum <- 0
    for (k in 1:length(partial_genelist)) {
      Ff_sum <- Ff_sum + partial_genelist[k] * (partial_genelist[k]-1)
    }
    Ff_sum <- as.numeric(Ff_sum)
    #Y <- ncol(data)
    N <- sum(partial_genelist)
    D <- N %% Y
    H <- N %% num

    CWrel <- (Y*(N-D+Ff_sum)-N^2+D^2)/(Y*(H^2+num*(N-H)-D)-N^2+D^2)
    return(CWrel)
  }

### End
### ------------------------------------------------------------------------ ###

