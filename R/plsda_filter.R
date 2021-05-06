options(warn = -1)
OPLSDA_test <- function(mat, lab, cutoff = 1) {
    X <- mat
    Y <- as.factor(lab)
    # !!! adding a FALSE for scale.
    # For the following PLS-DA model, the samples in row and variables in column.
    oplsda <- ropls::opls(X, Y, permI = 100, predI = 2, scaleC = "standard", fig.pdfC = FALSE)

    #res <- oplsda$vipVn
    res <- oplsda@vipVn
    cpds <- data.frame(CompoundName = names(res), VIP = res)
    cpds_filter <- (cpds$VIP > cutoff)
    return(cpds_filter)
}
