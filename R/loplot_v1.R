#' @title Loplot
#' @description Loplot, a visible figure of signal correction loplot
#' provide the visible figure of QC-based correction. See the details at the following References.
#' @param x the file before QC-RLS correction.
#' @param z the file after QC-RLS correction.
#' @param i a index for the name of variable.
#' @usage loplot(x,z,i)
#' @export
#' @references statTarget: a streamlined tool for signal drift correction
#' and interpretations of quantitative mass spectrometry-based omics data.
#' Luan H, Ji F, Chen Y, Cai Z. 2018, Analytica Chimica Acta.
loplot <- function(x,z,i){
  # x is the loess
  cn <- colnames(x)
  qcid <- grep("QC",cn)

  RSD30_CV=paste(rownames(x)[i],"_", i,".png", sep="")
  dirout.loplot <- paste(getwd(), "/statTarget/shiftCor/After_shiftCor/loplot", sep="")
  dir.create(dirout.loplot)

  png(paste(dirout.loplot,RSD30_CV,sep="/"))
  graphics::layout(matrix(1:2,nrow=2))

  numY <- 1:dim(x)[2]
  graphics::plot(numY,x[i,],pch=19,col="yellow",ylab = c("Intensity"),
                 xlab = c("Injection Order"), main = "Raw Peak")
  points(qcid,x[i,qcid],pch=19,col="blue")
  legend("top", c("Sample", "QC"),col=c('yellow', 'blue'),
         lty=1,pch= 19,bty="n", cex=0.75,horiz = TRUE)
  #lines(qcid,x[i,qcid],col=rgb(0,0,0,0.3),lwd=4)
  loe <- loess(x[i,qcid]~qcid)
  points(numY,predict(loe,numY),type='l',col=rgb(0,0,0,0.3),lwd=4)

  graphics::plot(numY,z[i,],pch=19,col="yellow",ylab = c("Intensity"),
                 xlab = c("Injection Order"),main = "Corrected Peak")
  points(qcid,z[i,qcid],pch=19,col="blue")
  #abline(h = 1, type='l',col=rgb(0,0,0,0.3),lwd=4)
  #lines(qcid,z[i,qcid],col=rgb(0,0,0,0.3),lwd=4)
  #loe_n <- loess(z[i,qcid]~qcid)
  #points(numY,predict(loe_n,numY),type='l',col=rgb(0,0,0,0.3),lwd=4)
  legend("top", c("Sample", "QC"),col=c('yellow', 'blue'),
         lty=1,pch= 19,bty="n", cex=0.75,horiz = TRUE)
  dev.off()
}
