#' @importFrom ropls opls
options(warn = -1)
OPLSDA_C <- function(mat, lab) {
  X <- mat
  Y <- as.factor(lab)

  # !!! adding a FALSE for scale.
  # For the following PLS-DA model, the samples in row and variables in column.
  #oplsda <- opls(X, Y, permI = 100, predI = 2, scaleC = "standard", crossvalI=2, plot = FALSE)
  oplsda <- ropls::opls(X, Y, permI = 100, predI = 2, scaleC = "standard", crossvalI=2, fig.pdfC = FALSE)
  #res <- oplsda$vipVn
  res <- oplsda@vipVn
  cpds <- data.frame(CompoundName = names(res), VIP = res)

  cpds_filter <- cpds[order(as.numeric(cpds[,2]), decreasing=T),1]

  return(cpds_filter)
}
##########################################################
