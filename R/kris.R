#' KRIS: Keen and Reliable Interface Subroutines for Bioinformatic Analysis
#'
#' Tools for common steps in the bioinformatic analysis of numeric and SNP
#' data.
#'
#' Functions:
#' \itemize{
#' \item \code{\link{cal.pc.linear}}: linear principal component analysis
#' (PCA).
#' \item \code{\link{cal.pc.projection}}: PCA of controls, with cases projected
#' onto it.
#' \item \code{\link{fst.hudson}}: average fixation index (Fst) between two
#' groups.
#' \item \code{\link{fst.each.snp.hudson}}: Fst for each SNP between two
#' groups.
#' \item \code{\link{plot3views}}: scatter plots of 3 variables in three views.
#' \item \code{\link{read.bed}} and \code{\link{write.bed}}: read and write
#' binary PLINK files.
#' \item \code{\link{rubikclust}}: unsupervised clustering to find rough
#' structure and outliers.
#' \item \code{\link{xxt}}: \code{X \%*\% t(X)} for large matrices, computed in
#' blocks.
#' }
#'
#' Example data: \code{\link{simsnp}} and \code{\link{sample_labels}}, loaded
#' with \code{data(example_SNP)}.
#'
#' KRIS uses \pkg{rARPACK} for the eigen decomposition.
#' @keywords internal
#' @import rARPACK
"_PACKAGE"
#> [1] "_PACKAGE"

NULL
