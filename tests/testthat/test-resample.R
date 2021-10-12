test_that("bootsample", {
  n <- 1000
  mu.X <- c(1, 2, 3)
  Sigma.X <- matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
  X.complete.cont <- MASS::mvrnorm(n, mu.X, Sigma.X)
  ls_boot <- bootsample(X.complete.cont, 4)
  expect_equal(length(ls_boot), 4)
  expect_equal(nrow(ls_boot[[1]]), n)
  expect_equal(ncol(ls_boot[[1]]), 3)
})

test_that("jacksample", {
  n <- 1000
  mu.X <- c(1, 2, 3)
  Sigma.X <- matrix(c(9, 3, 2, 3, 4, 0, 2, 0, 1), nrow = 3)
  X.complete.cont <- MASS::mvrnorm(n, mu.X, Sigma.X)
  ls_jack <- jacksample(X.complete.cont, 4)
  expect_equal(length(ls_jack), 4)
  expect_equal(nrow(ls_jack[[1]]), 3 * n / 4)
  expect_equal(ncol(ls_jack[[1]]), 3)
})
