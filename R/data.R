#' Simulated SNP dataset
#'
#' 2,000 independent SNPs for 103 simulated individuals: two populations of 50
#' individuals each, plus 3 outliers. The genetic distance between the
#' populations was set to Fst = 0.01, as in Balding and Nichols (1995).
#'
#' \describe{
#'   \item{snp}{A numeric matrix with 103 rows (individuals) and 2,000 columns
#'   (SNPs), in additive coding (0, 1, 2).}
#'   \item{snp.info}{A data.frame with one row per SNP and the columns chr
#'   (chromosome), ID (SNP ID), GD (genetic distance), position, allele1 and
#'   allele2.}
#'   \item{ind.info}{A data.frame with one row per individual and the columns
#'   FamID (family ID), IndID (individual ID), PatID (father ID), MatID
#'   (mother ID), sex and phenotype.}
#' }
#'
#' @name simsnp
#' @docType data
#' @format A list with 3 elements: \code{snp}, \code{snp.info} and
#' \code{ind.info}.
#' @seealso \code{\link{sample_labels}}
#' @usage data(example_SNP)
#' @keywords simsnp
#' @references
#' Balding, D.J., and Nichols, R.A. (1995). A method for quantifying
#' differentiation between populations at multi-allelic loci and its
#' implications for investigating identity and paternity. Genetica 96, 3-12.
"simsnp"

#' Population labels for simsnp
#'
#' The population label of each of the 103 individuals in
#' \code{\link{simsnp}}: "pop1", "pop2" or "outlier3".
#'
#' @name sample_labels
#' @docType data
#' @format A vector with 103 elements.
#' @seealso \code{\link{simsnp}}
#' @usage data(example_SNP)
#' @keywords sample_labels
"sample_labels"

NULL
