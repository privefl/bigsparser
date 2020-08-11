################################################################################

test_that("$add_columns() works with a dgCMatrix", {

  library(Matrix)
  spmat <- rsparsematrix(30, 30, density = 0.1, symmetric = TRUE)

  library(bigsparser)
  sfbm  <- as_SFBM(spmat[, 1:10])
  sfbm$add_columns(spmat[, 11:20], offset_i = 0)
  sfbm$add_columns(spmat[, 21:30], offset_i = 0)
  expect_equal(dim(sfbm), c(30, 30))
  expect_equal(sfbm$p, as(spmat, "dgCMatrix")@p)

  b <- runif(ncol(sfbm))
  expect_equal(sp_prodVec(sfbm, b), as.vector(spmat %*% b))
  expect_equal(sp_solve_sym(sfbm, b, add_to_diag = 1e-4),
               as.vector(solve(spmat + Diagonal(ncol(spmat), 1e-4), b)))
})

################################################################################

test_that("$add_columns() works with a dsCMatrix", {

  library(Matrix)
  spmat <- rsparsematrix(30, 30, density = 0.1, symmetric = TRUE)

  library(bigsparser)
  sfbm  <- as_SFBM(spmat)
  sfbm$add_columns(spmat, offset_i = 0)
  expect_equal(dim(sfbm), c(30, 60))
  sfbm2 <- as_SFBM(spmat2 <- cbind(spmat, spmat))
  expect_equal(sfbm$p, sfbm2$p)
  expect_identical(readBin(sfbm$sbk,  what = 1, n = 1000),
                   readBin(sfbm2$sbk, what = 1, n = 1000))

  sfbm$add_columns(spmat, offset_i = nrow(spmat))
  expect_equal(dim(sfbm), c(60, 90))
  sfbm3 <- as_SFBM(Matrix::bdiag(spmat2, spmat))
  expect_equal(sfbm$p, sfbm3$p)
  expect_identical(readBin(sfbm$sbk,  what = 1, n = 1000),
                   readBin(sfbm3$sbk, what = 1, n = 1000))
  # rbind(matrix(readBin(sfbm$sbk,  what = 1, n = 1000), nrow = 2),
  #       matrix(readBin(sfbm3$sbk, what = 1, n = 1000), nrow = 2))
})

################################################################################
