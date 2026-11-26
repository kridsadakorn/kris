#' (Internal) Replace missing values
#'
#' Used internally by the PCA functions.
#'
#' @param X A vector.
#' @param missing The value that marks a missing entry. Default NA.
#' @param rep A vector of replacement values, the same length as X.
#'
#' @return X with its missing entries replaced by the matching values of rep.
#'
replace.missing <- function(X,missing=NA,rep){
  if (is.na(missing)){
    idx = which(is.na(X))
    X[idx] = rep[idx]
  }else{
    idx = which(X == missing)
    X[idx] = rep[idx]
  }
  return(X)
}
