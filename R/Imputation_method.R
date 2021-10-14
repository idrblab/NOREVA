

### ---------------------------- Method - 03 --------------------------------###
## Perform K-nearest neighbor imputation (knn)

knn<-function(x,k){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
result<-impute.knn(filterdata, k, rowmax = 0.5, colmax = 0.8, maxp = 1500)
 cObs <- result$data
 return(cObs) 
}


### ---------------------------- Method - 05 --------------------------------###
## Perform Zero imputation (zero)

zero<-function(x){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata[is.na(filterdata)]<-0
  cObs <- filterdata
  return(cObs) 
}

### ---------------------------- Method - 06 --------------------------------###
## Perform Background imputation (back)

back<-function(x){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata<-x
  #filterdata[is.na(filterdata)] <- min(filterdata,na.rm=TRUE)
  filterdata[is.na(filterdata)] <- min(filterdata[filterdata>0],na.rm=TRUE)/2
  cObs <- filterdata
  return(cObs) 
}




### ---------------------------- Method - 08 --------------------------------###
## Column Median Imputation


ImputMedian <- function(x){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata<-x
  cObs <- rbind(NULL,NULL)
  i <- 3
  for (i in 1:dim(filterdata)[1]){
    filterdata[i,][is.na(filterdata[i,])] <- median(filterdata[i,],na.rm=TRUE)
    
  }
  cObs <- filterdata
  return(cObs)
}




### ---------------------------- Method - 09 --------------------------------###
## Column Minimum Imputation


ImputMean <- function(x){
  filterdata<-t(x)
  x<-filterdata
  x<-x[apply(x, 1, function(y) !all(is.na(y))),]
  filterdata<-x
  cObs <- rbind(NULL,NULL)
  i <- 3
  for (i in 1:dim(filterdata)[1]){
    filterdata[i,][is.na(filterdata[i,])] <- mean(filterdata[i,],na.rm=TRUE)
    
  }
  cObs <- filterdata
  return(cObs)
}















