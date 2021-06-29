################################################################################

test_that("can use SFBMs from old versions", {

  spmat <- Matrix::sparseMatrix(1, 2, x = 3, dims = c(3, 3))
  # saveRDS(bigsparser::as_SFBM(spmat), "inst/testdata/old_sfbm.rds", version = 2)
  test_file <- system.file("testdata/old_sfbm.rds", package = "bigsparser")
  if (file.exists(test_file)) {
    X <- readRDS(test_file)
    expect_false(identical(X$address, methods::new("externalptr")))
    expect_equal(sp_prodVec(X, rep(1, 3)), c(3, 0, 0))
  }
})

################################################################################

test_that("internal function work properly", {

  corr2 <- Matrix::Diagonal(4, 0:3)
  corr2[4, 2] <- 5
  corr2[1, 4] <- 6
  corr2[3, 4] <- 7

  test1 <- bigsparser:::range_col(corr2@p, corr2@i)
  expect_equal(test1, list(c(-1, 1, 2, 0), c(-2, 3, 2, 3)))

  test2 <- bigsparser:::write_val_compact(
    tmp <- tempfile(), corr2@p, corr2@i, corr2@x,
    offset_p = 0, offset_i = 0, symmetric = FALSE)
  expect_equal(test2, c(test1[1], list(c(0, 0, 3, 4, 8))))

  expect_identical(readBin(tmp, what = 1, n = 100),
                   c(1, 0, 5, 2, 6, 0, 7, 3))


  (corr3 <- Matrix::forceSymmetric(corr2))

  test3 <- bigsparser:::range_col_sym(corr3@p, corr3@i)
  expect_equal(test3, list(c(3, 1, 2, 0), c(3, 1, 3, 3)))

  test4 <- bigsparser:::write_val_compact(
    tmp2 <- tempfile(), corr3@p, corr3@i, corr3@x,
    offset_p = 0, offset_i = 0, symmetric = TRUE)
  expect_equal(test4, c(test3[1], list(c(0, 1, 2, 4, 8))))

  expect_identical(readBin(tmp2, what = 1, n = 100),
                   c(6, 1, 2, 7, 6, 0, 7, 3))

  COL_RANGE <- function(x) {
    ind <- which(x != 0)
    `if`(length(ind) > 0, range(ind) - 1L, c(-1, -2))
  }
  spmat <- Matrix::rsparsematrix(1000, 1000, 0.001)
  test5 <- bigsparser:::range_col(spmat@p, spmat@i)
  expect_equal(do.call("rbind", test5), apply(spmat, 2, COL_RANGE))

  spmat2 <- Matrix::rsparsematrix(1000, 1000, 0.001, symmetric = TRUE)
  test6 <- bigsparser:::range_col_sym(spmat2@p, spmat2@i)
  expect_equal(do.call("rbind", test6), apply(spmat2, 2, COL_RANGE))
})

################################################################################

test_that("can create an SFBM_compact from a dgCMatrix", {

  spmat <- Matrix::rsparsematrix(1000, 1000, 0.01)
  X <- as_SFBM(spmat, compact = TRUE)
  expect_true(all(X$p >= spmat@p))
  expect_equal(dim(X), spmat@Dim)

  expect_false(X$is_saved)
  X <- X$save()
  expect_true(X$is_saved)
  X2 <- readRDS(X$rds)
  expect_equal(X2$p, X$p)
  expect_equal(X2$sbk, X$sbk)
})

################################################################################

test_that("can create an SFBM_compact from a dsCMatrix", {

  spmat0 <- Matrix::rsparsematrix(1000, 1000, 0.01, symmetric = TRUE)
  X <- as_SFBM(spmat0, compact = TRUE)
  spmat <- as(spmat0, "dgCMatrix")
  X2 <- as_SFBM(spmat, compact = TRUE)
  expect_identical(readBin(X$sbk,  what = 1, n = 1e6),
                   readBin(X2$sbk, what = 1, n = 1e6))
  expect_equal(dim(X), spmat@Dim)
})

################################################################################

test_that("products between an SFBM_compact and a vector work", {

  spmat <- Matrix::rsparsematrix(1000, 1000, 0.1)
  X <- as_SFBM(spmat, compact = TRUE)
  y <- runif(1000)
  expect_equal(sp_prodVec(X, y), as.vector(spmat %*% y))
  expect_equal(sp_cprodVec(X, y), as.vector(Matrix::crossprod(spmat, y)))
})

################################################################################
