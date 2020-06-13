################################################################################

test_that("$add_columns() works", {

  library(Matrix)
  spmat <- rsparsematrix(30, 30, density = 0.1, symmetric = TRUE)

  library(bigsparser)
  sfbm  <- as_SFBM(spmat[, 1:10])
  sfbm$add_columns(spmat[, 11:20], offset_i = 0)
  sfbm$add_columns(spmat[, 21:30], offset_i = 0)
  expect_equal(dim(sfbm), c(30, 30))

  b <- rep(1, ncol(sfbm))
  expect_equal(sp_prodVec(sfbm, b), as.vector(spmat %*% b))
  expect_equal(sp_solve_sym(sfbm, b, add_to_diag = 1e-4),
               as.vector(solve(spmat + Diagonal(ncol(spmat), 1e-4), b)))
})

################################################################################
