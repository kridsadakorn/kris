#' Multiply a large matrix by its transpose
#'
#' @description Computes \code{X \%*\% t(X)} by splitting the columns of X into
#' blocks of \code{window.size} columns and adding up the products of the
#' blocks. This "divide and conquer" approach speeds up the computation for
#' matrices with many columns.
#'
#' @param X A numeric matrix.
#' @param window.size Number of columns in each block. Must be a positive
#' integer no larger than \code{ncol(X)}. Default 5.
#'
#' @return The matrix \code{X \%*\% t(X)}. If the input is invalid, the
#' function prints a message and returns NULL.
#'
#' @export
#'
#' @examples
#'
#'
#' #Create a random matrix
#' X <-matrix(runif(100), ncol=20)
#' R1  <- xxt(X)
#'
#' #Show the result (R1)
#' print(R1)
#' R2 <- X %*% t(X)
#'
#' #Show the result (R2)
#' print(R2)
#'
#'
#'

xxt <- function(X, window.size = 5){
  n.col = dim(X)[2]
  window.size = as.integer(window.size)
  if (window.size < 1){
    cat(paste0("In matrices.split.by.col(), window.size must be a positive integer\n"))
    return(NULL)
  }
  if (is.null(n.col)){
    cat(paste0("In matrices.split.by.col(), n.col is null.\n"))
    return(NULL)
  }
  if (n.col < window.size){
    cat(paste0("In matrices.split.by.col(), n.col is less than window.size\n"))
    return(NULL)
  }

  no.split = n.col %/% window.size

  split.matrices = lapply(1:no.split, function(col){
    X[,((col-1)*window.size+1):(col*window.size), drop=FALSE]})

  if ((no.split * window.size) < n.col){
    split.matrices[[length(split.matrices)+1]] = X[,(no.split*window.size+1):n.col, drop=FALSE]
  }

  #Can be parallel

  partial.mul = NULL
  for (i in 1:length(split.matrices)){
    partial.mul[[i]] = split.matrices[[i]] %*% t(split.matrices[[i]])
  }

  mul = Reduce(`+`, partial.mul)

  return(mul)
}





