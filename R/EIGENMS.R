#' @title EigenMS normalization
#' @description EigenMS is an adaptation of surrogate variable analysis, which identifies trends
#' attributable to bias by utilizing singular value decomposition on model residuals.
#' See the details at the following References.
#' @param data Input matrix of data
#' @param label Input the label of data
#' @return A structure with multiple components
#' \describe{
#'   \item{m_logInts}{number of metabolites x number of samples
#'            matrix of expression data with no missing values}
#'   \item{grps}{either a single factor indicating the treatment group of
#'           each sample i.e. [1 1 1 1 2 2 2 2...]
#'           or a data frame of factors}
#'   \item{m_ints_eig1}{First portion of EigenMS: Identify eigentrends attributable to bias}
#'   \item{m_ints_norm1}{Eliminate the effects of systematic bias identified in eig_norm1()}
#'   \item{mm_eigenMS}{matrix of normalized abundances, no extra columns}
#'}
#' @importFrom ProteoMM eig_norm1
#' @importFrom ProteoMM eig_norm2
#' @usage EIGENMS(data,label)
#' @examples
#' \donttest{data(mm_metabolites)}
#' \donttest{head(mm_metabolites)}
#' # different from parameter names as R uses outer name spaces
#' # if variable is undefined
#' \donttest{intsCols = 8:13}
#' \donttest{metaCols = 1:7} # reusing this variable
#' \donttest{m_logInts = make_intencities(mm_metabolites, intsCols)}  # will reuse the name
#' \donttest{m_logInts = convert_log2(m_logInts)}
#' # 3 samples for CG and 3 for mCG
#' \donttest{grps = as.factor(c('CG','CG','CG', 'mCG','mCG','mCG'))}
#' \donttest{mm_eigenMS = EIGENMS(m_logInts,grps)}
#' @references 1. Metabolomics data normalization with EigenMS.
#' Karpievitch YK, Nikolic SB, Wilson R, Sharman JE, Edwards LM. 2014, PLoS ONE.
#' @references 2. Normalization of peak intensities in bottom-up MS-based proteomics
#' using singular value decomposition.
#' Karpievitch YV, Taverner T et al. 2009, Bioinformatics.

#' @export

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
