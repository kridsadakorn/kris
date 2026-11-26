test_that("ISSUE-006: plot3views restores par and accepts matrix input", {
  pdf(file.path(tempdir(), "issue006.pdf"))
  on.exit(dev.off())
  before <- par("mfrow", "mar")
  X <- matrix(seq_len(30), ncol = 3)
  labels <- matrix(rep(c("a", "b"), 5), ncol = 1)
  expect_null(plot3views(X, labels))
  expect_equal(par("mfrow", "mar"), before)
})
