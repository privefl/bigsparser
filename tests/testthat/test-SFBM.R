################################################################################

test_that("can create an SFBM from a dgCMatrix", {

  spmat <- Matrix::rsparsematrix(1000, 1000, 0.01)
  X <- as_SFBM(spmat)
  expect_identical(X$p, spmat@p)
  expect_equal(dim(X), spmat@Dim)
  expect_equal(X$nval, length(spmat@x))

  expect_equal(file.size(X$sbk), length(spmat@x) * 16)
  con <- file(X$sbk, open = "rb")
  i <- rep(0, length(spmat@x))
  x <- rep(0, length(spmat@x))
  for (k in seq_along(x)) {
    i[k] <- readBin(con, n = 1, what = double())
    x[k] <- readBin(con, n = 1, what = double())
  }
  expect_equal(i, spmat@i)
  expect_identical(x, spmat@x)
  close(con)

  expect_false(X$is_saved)
  X <- X$save()
  expect_true(X$is_saved)
  X2 <- readRDS(X$rds)
  expect_equal(X2$p, X$p)
  expect_equal(X2$sbk, X$sbk)
})

################################################################################

test_that("can create an SFBM from a dsCMatrix", {

  spmat0 <- Matrix::rsparsematrix(1000, 1000, 0.01, symmetric = TRUE)
  X <- as_SFBM(spmat0)
  spmat <- as(spmat0, "dgCMatrix")
  X2 <- as_SFBM(spmat)
  expect_identical(readBin(X$sbk,  what = 1, n = 1e6),
                   readBin(X2$sbk, what = 1, n = 1e6))
  expect_equal(X$p, spmat@p)
  expect_equal(dim(X), spmat@Dim)
  expect_equal(X$nval, length(spmat@x))
  expect_equal(file.size(X$sbk), length(spmat@x) * 16)

  expect_false(X$is_saved)
  X <- X$save()
  expect_true(X$is_saved)
  X2 <- readRDS(X$rds)
  expect_equal(X2$p, X$p)
  expect_equal(X2$sbk, X$sbk)
})

################################################################################

test_that("products between an SFBM and a vector work", {

  spmat <- Matrix::rsparsematrix(1000, 1000, 0.1)
  X <- as_SFBM(spmat, compact = sample(c(TRUE, FALSE), 1))
  y <- runif(1000)
  expect_equal(sp_prodVec(X, y), as.vector(spmat %*% y))
  expect_equal(sp_cprodVec(X, y), as.vector(Matrix::crossprod(spmat, y)))
})

################################################################################

test_that("SFBM_corr_compact works", {

 replicate(5, {

   spmat <- Matrix::rsparsematrix(1000, 1000, 0.1)
   diag(spmat) <- max(abs(spmat))
   corr <- Matrix::cov2cor(Matrix::forceSymmetric(spmat))
   X2 <- as_SFBM_corr_compact(corr)

   y <- runif(1000)
   expect_equal(sp_prodVec(X2, y), as.vector(corr %*% y), tolerance = 1e-4)
   expect_equal(sp_cprodVec(X2, y), as.vector(Matrix::crossprod(corr, y)),
                tolerance = 1e-4)

   I <- sample(1000, 500); J <- sample(1000, 500)
   expect_lt(max(abs(X2[I, J] - corr[I, J])), 1e-4)
   expect_null(bigassertr::assert_int(X2[I, J]@x * 32767))
 })
})

################################################################################
