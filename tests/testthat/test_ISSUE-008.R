test_that("ISSUE-008: read.bed returns what write.bed wrote, for any size", {
  set.seed(8)
  for (sz in list(c(10, 6), c(5, 3), c(7, 1), c(4, 4), c(9, 13), c(1, 25))) {
    n <- sz[1]; p <- sz[2]
    snp <- matrix(sample(c(0, 1, 2, NA), n * p, TRUE), n, p)
    obj <- list(snp = snp,
                snp.info = data.frame(1, paste0("m", 1:p), 0, 1:p, "A", "T"),
                ind.info = data.frame(1:n, 1:n, 0, 0, 1, 1))
    f <- file.path(tempdir(), "issue008")
    write.bed(obj, f)
    back <- read.bed(paste0(f, ".bed"), paste0(f, ".bim"), paste0(f, ".fam"))
    expect_equal(unname(back$snp), snp, info = paste(n, "x", p))
    file.remove(paste0(f, c(".bed", ".bim", ".fam")))
  }
})

test_that("ISSUE-008: read.bed matches the example data", {
  bed <- system.file("extdata", "example_SNP.bed", package = "KRIS")
  bim <- system.file("extdata", "example_SNP.bim", package = "KRIS")
  fam <- system.file("extdata", "example_SNP.fam", package = "KRIS")
  snp <- read.bed(bed, bim, fam)
  data(example_SNP)
  expect_equal(dim(snp$snp), c(103, 2000))
  expect_equal(unname(snp$snp), unname(simsnp$snp))
})

test_that("ISSUE-008: xxt works on a one-row matrix", {
  X <- matrix(1:10, nrow = 1)
  expect_equal(xxt(X), X %*% t(X))
})

test_that("ISSUE-008: cal.pc.linear gives the same result with XXT TRUE or FALSE", {
  data(example_SNP)
  a <- cal.pc.linear(simsnp$snp, no.pc = 3, XXT = TRUE)
  b <- cal.pc.linear(simsnp$snp, no.pc = 3, XXT = FALSE)
  expect_equal(a$evalue, b$evalue, tolerance = 1e-6)
  expect_equal(abs(a$PC), abs(b$PC), tolerance = 1e-6)
})

test_that("ISSUE-008: no.pc = 1 works", {
  data(example_SNP)
  expect_equal(dim(cal.pc.linear(simsnp$snp, no.pc = 1)$PC), c(103, 1))
  st <- rep(c(1, 2), length.out = nrow(simsnp$snp))
  expect_equal(dim(cal.pc.projection(simsnp$snp, st, no.pc = 1)$PC), c(103, 1))
})
