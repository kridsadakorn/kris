context("Read PLINK format")

test_that("Read BED file",{
  library(kris)
  bed <- system.file("extdata", "example_SNP.bed", package="kris")
  bim <- system.file("extdata", "example_SNP.bim", package="kris")
  fam <- system.file("extdata", "example_SNP.fam", package="kris")
  snp <- read.bed(bed, bim, fam )

  expect_length(ls(snp), 3)
  expect_type(ls(snp), "character")
})
