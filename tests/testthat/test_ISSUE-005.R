test_that("ISSUE-005: xxt works for any valid window.size", {
  X <- matrix(seq_len(100) / 7, ncol = 20)
  for (w in c(1, 5, 7, 11, 20)) {
    expect_equal(xxt(X, window.size = w), X %*% t(X))
  }
})
