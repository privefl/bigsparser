################################################################################

test_that("$add_columns() works with a dgCMatrix", {

  spmat <- Matrix::rsparsematrix(30, 30, density = 0.1, symmetric = TRUE)

  sfbm  <- as_SFBM(spmat[, 1:10])
  sfbm$add_columns(spmat[, 11:20], offset_i = 0)
  sfbm$add_columns(spmat[, 21:30], offset_i = 0)
  expect_equal(dim(sfbm), c(30, 30))
  expect_equal(sfbm$p, as(spmat, "dgCMatrix")@p)

  b <- runif(ncol(sfbm))
  expect_equal(sp_prodVec(sfbm, b), as.vector(spmat %*% b))
  expect_equal(sp_solve_sym(sfbm, b, add_to_diag = 1e-4),
               as.vector(solve(spmat + Matrix::Diagonal(ncol(spmat), 1e-4), b)))
})

################################################################################

test_that("$add_columns() works with a dsCMatrix", {

  spmat <- Matrix::rsparsematrix(30, 30, density = 0.1, symmetric = TRUE)

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

test_that("$add_columns() (compact) works with a dgCMatrix", {

  spmat <- Matrix::rsparsematrix(30, 30, density = 0.1, symmetric = TRUE)

  sfbm  <- as_SFBM(spmat[, 1:10], compact = TRUE)
  sfbm$add_columns(spmat[, 11:20], offset_i = 0)
  sfbm$add_columns(spmat[, 21:30], offset_i = 0)
  expect_equal(dim(sfbm), c(30, 30))
  sfbm2 <- as_SFBM(spmat, compact = TRUE)
  expect_identical(sfbm$p, sfbm2$p)
  expect_identical(sfbm$first_i, sfbm2$first_i)
  expect_identical(readBin(sfbm$sbk, 1, 1000), readBin(sfbm2$sbk, 1, 1000))

  b <- runif(ncol(sfbm))
  expect_equal(sp_prodVec(sfbm, b), as.vector(spmat %*% b))
  expect_equal(sp_solve_sym(sfbm, b, add_to_diag = 1e-4),
               as.vector(solve(spmat + Matrix::Diagonal(ncol(spmat), 1e-4), b)))
})

################################################################################

test_that("$add_columns() (compact) works with a dsCMatrix", {

  spmat <- Matrix::rsparsematrix(30, 30, density = 0.1, symmetric = TRUE)

  sfbm  <- as_SFBM(spmat, compact = TRUE)
  sfbm$add_columns(spmat, offset_i = 0)
  expect_equal(dim(sfbm), c(30, 60))
  spmat2 <- cbind(spmat, spmat)
  sfbm2 <- as_SFBM(spmat2, compact = TRUE)
  expect_equal(sfbm$p, sfbm2$p)
  expect_equal(sfbm$first_i, sfbm2$first_i)
  expect_identical(readBin(sfbm$sbk,  what = 1, n = 1000),
                   readBin(sfbm2$sbk, what = 1, n = 1000))

  sfbm$add_columns(spmat, offset_i = nrow(spmat))
  expect_equal(dim(sfbm), c(60, 90))
  spmat3 <- Matrix::bdiag(spmat2, spmat)
  sfbm3 <- as_SFBM(spmat3, compact = TRUE)
  expect_equal(sfbm$p, sfbm3$p)
  expect_equal(sfbm$first_i, sfbm3$first_i)
  expect_identical(readBin(sfbm$sbk,  what = 1, n = 1000),
                   readBin(sfbm3$sbk, what = 1, n = 1000))

  b <- runif(ncol(sfbm3))
  expect_equal(sp_prodVec(sfbm3, b), as.vector(spmat3 %*% b))
})

################################################################################

test_that("$add_columns() works with large data", {

  skip_on_cran()
  skip_on_covr()

  spmat <- Matrix::rsparsematrix(20e3, 20e3, nnz = 1e7, symmetric = TRUE)
  sfbm  <- as_SFBM(spmat)

  beg <- readBin(sfbm$sbk, what = 0, n = 1e6)

  all_time <- replicate(15, {
    # cat(".")
    time <- system.time(sfbm$add_columns(spmat, nrow(sfbm)))[3]
    # print(file.size(sfbm$backingfile))
    time
  })
  # plot(all_time)
  expect_equal(dim(sfbm), 16 * rep(20e3, 2))
  expect_equal(sfbm$nval, 16 * length(as(spmat, "dgCMatrix")@x))
  expect_identical(readBin(sfbm$sbk, what = 0, n = 1e6), beg)
})

################################################################################
