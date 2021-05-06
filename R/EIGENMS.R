#' title EigenMS normalization
#' description EigenMS estimates and preserves fixed effects.
#' Random effects may be attempted in the future.
#' @param data This is the first ones argument
#' @param label This is the second ones argument
#' @importFrom ProteoMM eig_norm1
#' @importFrom ProteoMM eig_norm2

EIGENMS <- function(data, label) {

  a <- data
  ddata<-a
  m_logInts = ddata
  grps = as.factor(label)
  m_prot.info = cbind(rownames(ddata),rownames(ddata))
  m_ints_eig1 = suppressWarnings(suppressMessages(eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)))
  m_ints_norm1 = suppressWarnings(suppressMessages(eig_norm2(rv=m_ints_eig1)))
  eigenMS <-m_ints_norm1$norm_m
  return(eigenMS)

}
