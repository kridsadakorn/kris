#' Average fixation index (Fst) between two groups
#'
#' @description Computes the average Fst over all SNPs between two groups of
#' individuals with the Hudson method (Hudson et al. 1992). The average is
#' taken as a ratio of averages, as recommended by Bhatia et al. (2013).
#'
#' @param X A matrix of SNPs in additive coding (0, 1, 2), with individuals in
#' rows and SNPs in columns.
#' @param idx.p1 Row indices of the individuals in the first group.
#' @param idx.p2 Row indices of the individuals in the second group.
#'
#' @return A single Fst value.
#'
#' @references
#' Bhatia, G., Patterson, N., Sankararaman, S., and Price, A.L. (2013).
#' Estimating and interpreting FST: The impact of rare variants. Genome Res. 23,
#' 1514-1521.
#'
#' Hudson, R.R., Slatkin, M., and Maddison, W.P. (1992). Estimation of levels of
#' gene flow from DNA sequence data. Genetics 132, 583-589.
#'
#' @export
#'
#' @seealso \code{\link{fst.each.snp.hudson}}
#'
#' @examples
#'
#'
#' #Load simulated dataset
#' data(example_SNP)
#'
#' idx1 <- which(sample_labels == 'pop1')
#' idx2 <- which(sample_labels == 'pop2')
#' fst <- fst.hudson(simsnp$snp, idx1, idx2)
#'
#' #Print out the Fst value between 'pop1' and 'pop2'
#' print(fst)
#'

fst.hudson <-function(X, idx.p1, idx.p2){
  prestep.fst.one.marker <- function(alleles,idx.p1,idx.p2){

    #Pop 1
    G = alleles[idx.p1]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n1=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n1
    #p.B=(no.BB*2 + no.AB)/n1
    #p1 = min(p.A,p.B,na.rm = T)
    p1 = p.A

    #Pop 2
    G = alleles[idx.p2]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n2=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n2
    #p.B=(no.BB*2 + no.AB)/n2
    #p2 = min(p.A,p.B,na.rm = T)
    p2 = p.A


    N = (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    D = p1*(1-p2) + p2*(1-p1)

    return(c(N,D))
  }

  set.fst = apply(X,2,prestep.fst.one.marker,idx.p1=idx.p1, idx.p2=idx.p2)

  # Bhatia, et al 2013, "This is the basis of our recommendation that FST be estimated as a ratio of averages."

  fst = mean(set.fst[1,],na.rm=T) / mean(set.fst[2,],na.rm=T)
  #fst = mean(set.fst[1,]/set.fst[2,],na.rm = T)

  return(fst)
}

#' Fixation index (Fst) for each SNP between two groups
#'
#' @description Computes Fst for every SNP between two groups of individuals
#' with the Hudson method (Hudson et al. 1992; Bhatia et al. 2013).
#'
#' @param X A matrix of SNPs in additive coding (0, 1, 2), with individuals in
#' rows and SNPs in columns.
#' @param idx.p1 Row indices of the individuals in the first group.
#' @param idx.p2 Row indices of the individuals in the second group.
#'
#' @return A numeric vector with one Fst value per SNP (column of X).
#'
#' @references
#' Bhatia, G., Patterson, N., Sankararaman, S., and Price, A.L. (2013).
#' Estimating and interpreting FST: The impact of rare variants. Genome Res. 23,
#' 1514-1521.
#'
#' Hudson, R.R., Slatkin, M., and Maddison, W.P. (1992). Estimation of levels of
#' gene flow from DNA sequence data. Genetics 132, 583-589.
#'
#' @export
#'
#' @seealso \code{\link{fst.hudson}}
#'
#' @examples
#'
#'
#' #Load simulated dataset
#' data(example_SNP)
#'
#' idx1 <- which(sample_labels == 'pop1')
#' idx2 <- which(sample_labels == 'pop2')
#' fst.pairwise <- fst.each.snp.hudson(simsnp$snp, idx1, idx2)
#'
#' #Print out the Fst values of the first three SNPs between 'pop1' and 'pop2'
#' print(fst.pairwise[1:3])
#'

fst.each.snp.hudson <-function(X, idx.p1, idx.p2){
  prestep.fst.one.marker <- function(alleles,idx.p1,idx.p2){

    #Pop 1
    G = alleles[idx.p1]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n1=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n1
    #p.B=(no.BB*2 + no.AB)/n1
    #p1 = min(p.A,p.B,na.rm = T)
    p1 = p.A

    #Pop 2
    G = alleles[idx.p2]
    no.AA=length(which(G==0))
    no.AB=length(which(G==1))
    no.BB=length(which(G==2))
    n2=(no.AA+no.AB+no.BB)*2
    p.A=(no.AA*2 + no.AB)/n2
    #p.B=(no.BB*2 + no.AB)/n2
    #p2 = min(p.A,p.B,na.rm = T)
    p2 = p.A


    N = (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    D = p1*(1-p2) + p2*(1-p1)

    return(c(N,D))
  }

  set.fst = apply(X,2,prestep.fst.one.marker,idx.p1=idx.p1, idx.p2=idx.p2)
  #fst = mean(set.fst[1,],na.rm=T) / mean(set.fst[2,],na.rm=T)
  vec.fst = set.fst[1,]/set.fst[2,]

  return(vec.fst)
}

