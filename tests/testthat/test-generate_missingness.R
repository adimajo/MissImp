test_that("generate_missingness", {
  n <- 1000
  mu.X <- c(1, 2, 3)
  Sigma.X <- matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
  X.complete.cont <- MASS::mvrnorm(n, mu.X, Sigma.X)
  rs <- generate_miss(X.complete.cont, 0.5, mechanism = "MNAR2")
  expect_equal(length(rs), 3)
  expect_equal(nrow(rs$X.incomp), n)
  expect_equal(ncol(rs$X.incomp), 3)
})
