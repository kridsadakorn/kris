#' Linear principal component analysis (PCA)
#'
#' @description Computes principal components (PCs) from a numeric matrix or
#' from SNP data in additive coding (0, 1, 2). Missing values are replaced by
#' the column median before the analysis.
#'
#' @param X A numeric matrix with samples in rows and features in columns.
#' @param PCscore If TRUE (default), return PC scores. If FALSE, return the
#' eigenvectors.
#' @param no.pc Number of PCs to compute. If NA (default), all PCs are
#' computed.
#' @param data.type Type of data in X: "linear" (default) for numeric data or
#' "snp" for SNPs in additive coding. It sets how each column is scaled.
#' @param XXT If TRUE (default), compute the PCs from \code{X \%*\% t(X)}; if
#' FALSE, from X directly. Both give the same PCs, but TRUE is faster when X
#' has many columns.
#'
#' @return A list with:
#' \itemize{
#' \item \code{$PC}: a matrix of PCs, with samples in rows and PCs in columns.
#' \item \code{$evalue}: a vector of eigenvalues.
#' }
#'
#' @export
#' @importFrom stats median sd
#'
#' @import rARPACK
#'
#' @seealso \code{\link{cal.pc.projection}}
#' @examples
#'
#' #Load simulated dataset
#' data(example_SNP)
#'
#' #Using default parameters
#' PCs <- cal.pc.linear(simsnp$snp)
#' summary(PCs)
#'
#' #Preview $PC
#' print(PCs$PC[1:5,1:3])
#'
#' #Preview $evalue
#' print(PCs$evalue[1:3])
#'
#' plot3views(PCs$PC[,1:3], sample_labels)
#'
#' #Calculate PCs without PC scores
#'
#' PCs <- cal.pc.linear(simsnp$snp, PCscore = FALSE)
#' summary(PCs)
#'
#' #Preview $PC
#' print(PCs$PC[1:5,1:3])
#'
#' #Preview $evalue
#' print(PCs$evalue[1:3])
#'
#' plot3views(PCs$PC[,1:3], sample_labels)
#'
#' #Calculate the top 3 PCs
#' PCs <- cal.pc.linear(simsnp$snp, no.pc = 3)
#' summary(PCs)
#'
#' #Preview $PC
#' print(PCs$PC[1:5,1:3])
#'
#' #Preview $evalue
#' print(PCs$evalue[1:3])
#'
#' plot3views(PCs$PC[,1:3], sample_labels)


cal.pc.linear <- function(X, PCscore = TRUE, no.pc = NA, data.type = "linear",
                          XXT = TRUE){

  #Resolve missing value by median
  #cat(paste0("Missing values are set as medians\n"))
  X.median = apply(X,2,median,na.rm=TRUE)
  missing.char=NA
  X = t(apply(X,1,replace.missing,missing=missing.char,rep=X.median))

  if (data.type == "snp"){
    XPi = (colSums(X)+1)/(2+(2*dim(X)[1]))
    XSD = (XPi * (1-XPi))^0.5
  }else{
    #Assume linear
    XSD = apply(X,2,sd)
  }

  XCM = colMeans(X)
  A = t((t(X)-XCM)/XSD)
  A = A[,colSums(is.na(A))<nrow(A)]

  if (is.na(no.pc)){
    no.pc = dim(A)[1]
  }

  if (XXT == TRUE){

    #To handle large matrix
    no.col = dim(A)[2]
    if (no.col <= 1500){
      AA=A %*% t(A)
    }else{
      AA = xxt(A,window.size = 1000)
    }

    UDV = svds(AA, k=no.pc)
    eigen.value = UDV$d
    PCs = UDV$u
  }else{
    UDV = svds(A, k=no.pc)
    # Singular values of A are the square roots of the eigenvalues of A.t(A)
    eigen.value = UDV$d^2
    PCs = UDV$u
  }


  if (PCscore == TRUE){
    #Calculate the PCscore and the scaled-up PCs
    B.T=diag(1/sqrt(eigen.value), nrow=length(eigen.value)) %*% t(PCs) %*% A
    PCs=A %*% t(B.T)
  }

  return(list("PC"=PCs,"evalue"=eigen.value))
}


#' PCA of controls with cases projected onto it
#'
#' @description Computes principal components (PCs) from the control group
#' only, then projects the case group onto those PCs. The disease status of
#' every individual is required. Missing values are replaced by the column
#' median before the analysis.
#'
#' @param X A numeric matrix with samples in rows and features in columns.
#' @param status A numeric vector with the disease status of each individual:
#' 1 for control, 2 for case. Individuals with any other value are left out of
#' the result.
#' @param individual_id Optional character vector of individual IDs.
#' @param labels Optional character vector of labels, such as population
#' names.
#' @param no.pc Number of PCs to compute. If NA (default), all PCs are
#' computed.
#' @param data.type Type of data in X: "linear" (default) or "snp".
#'
#' @return A list in which controls come first, then cases:
#' \itemize{
#' \item \code{$PC}: a matrix of PCs, with individuals in rows and PCs in
#' columns.
#' \item \code{$id}: the individual IDs, or NULL if \code{individual_id} was
#' not given.
#' \item \code{$label}: the labels, or NULL if \code{labels} was not given.
#' \item \code{$status}: "control" or "case" for each individual.
#' }
#'
#' @export
#' @import rARPACK
#'
#' @seealso \code{\link{cal.pc.linear}}
#'
#' @examples
#'
#' data(example_SNP)
#'
#' #Create a random list of disease status, 1 = Control and 2 = Case
#'
#' ind_status <- sample(c(1,2), size = length(sample_labels), replace = TRUE)
#'
#' PCs <- cal.pc.projection(simsnp$snp, status = ind_status,
#' labels = sample_labels)
#' summary(PCs)
#'
#' #Preview $PC
#' print(PCs$PC[1:5,1:3])
#'
#' #Preview $status
#' print(PCs$status[1:3])
#'
#' plot3views(PCs$PC[,1:3], PCs$label)
#'
#' #Calculate the top 3 PCs
#'
#' PCs <- cal.pc.projection(simsnp$snp, status = ind_status,
#' labels = sample_labels, no.pc = 3)
#' summary(PCs)
#'
#' #Preview $PC
#' print(PCs$PC[1:5,1:3])
#'
#' plot3views(PCs$PC[,1:3], PCs$label)
#'
#'
cal.pc.projection <- function(X, status, individual_id = NULL, labels = NULL,
                              no.pc = NA, data.type = "linear"){

  #Resolve missing value by median

  X.median = apply(X,2,median,na.rm=TRUE)
  missing.char=NA
  X = t(apply(X,1,replace.missing,missing=missing.char,rep=X.median))

  if (data.type == "snp"){
    XPi = (colSums(X)+1)/(2+(2*dim(X)[1]))
    XSD = (XPi * (1-XPi))^0.5
  }else{
    #Assume linear
    XSD = apply(X,2,sd)
  }

  X.con = X[which(status==1),]
  X.case = X[which(status==2),]

  if (!is.null(labels)){
    new.pop = labels[which(status==1)]
    tmp.pop = labels[which(status==2)]
    new.pop = c(new.pop,tmp.pop)
  }else{
    new.pop = NULL
  }

  if (!is.null(individual_id)){
    new.id = individual_id[which(status==1)]
    tmp.id = individual_id[which(status==2)]
    new.id = c(new.id,tmp.id)
  }else{
    new.id = NULL
  }

  X = X.con

  XPi = (colSums(X)+1)/(2+(2*dim(X)[1]))
  XSD = (XPi * (1-XPi))^0.5

  XCM = colMeans(X)
  A = t((t(X)-XCM)/XSD)
  A = A[,colSums(is.na(A))<nrow(A)]

  if (is.na(no.pc)){
    no.pc = dim(A)[1]
  }

  AA=A %*% t(A)
  evv = eigs_sym(AA, k=no.pc)
  eigen.value = evv$values
  PCs = evv$vectors
  B.T=diag(1/sqrt(eigen.value), nrow=length(eigen.value)) %*% t(PCs) %*% A

  PCcon=A %*% t(B.T)

  X = X.case

  XPi = (colSums(X)+1)/(2+(2*dim(X)[1]))
  XSD = (XPi * (1-XPi))^0.5

  XCM = colMeans(X)
  A = t((t(X)-XCM)/XSD)
  A = A[,colSums(is.na(A))<nrow(A)]

  PCcase=A %*% t(B.T)

  label.status=rep('control',dim(PCcon)[1])
  tmp.label=rep('case',dim(PCcase)[1])
  label.status = c(label.status,tmp.label)

  PC_project = rbind(PCcon,PCcase)

  return(list("PC"=PC_project,"id"=new.id,"label"=new.pop,"status"=label.status))
}

