################################################################################

test_that("can use SFBMs from old versions", {

  spmat <- Matrix::sparseMatrix(1, 2, x = 3, dims = c(3, 3))
  test_file <- test_path("testdata/old_sfbm.rds")
  X <- readRDS(test_file)
  X$backingfile <- sub("\\.rds$", ".sbk", test_file)
  expect_identical(X[], spmat)
})

################################################################################

test_that("sparse matrix accessors work", {

  spmat <- Matrix::Diagonal(5, c(1, 0:3))
  spmat[4, 3] <- 5
  spmat[1, 4] <- 6
  spmat[3, 4] <- 7

  list_ind_row <- list_ind_col <-
    list(2, 1:5, c(TRUE, FALSE, TRUE, FALSE, TRUE), -3, -(3:5), sample(5, replace = TRUE))

  for (compact in c(TRUE, FALSE)) {
    X <- as_SFBM(spmat, compact = compact)
    for (ind_row in list_ind_row) {
      for (ind_col in list_ind_col) {
        expect_identical(X[ind_row, ind_col],
                         spmat[ind_row, ind_col, drop = FALSE])
        expect_error(X[ind_row, ind_col, drop = FALSE], "not implemented")
      }
    }
    expect_identical(X[],  spmat)
    expect_identical(X[,], spmat)
  }

  spmat2 <- Matrix::rsparsematrix(5, 5, density = 0.5, symmetric = TRUE)
  diag(spmat2) <- max(abs(spmat2)); spmat3 <- Matrix::cov2cor(spmat2)
  X2 <- as_SFBM_corr_compact(spmat3)
  for (ind_row in list_ind_row) {
    for (ind_col in list_ind_col) {
      expect_lt(max(abs(X2[ind_row, ind_col] - spmat3[ind_row, ind_col, drop = FALSE])),
                1 / 32767)
      expect_error(X[ind_row, ind_col, drop = FALSE], "not implemented")
    }
  }
})

################################################################################

test_that("dense matrix accessors work", {

  spmat <- Matrix::Diagonal(5, c(1, 0:3))
  spmat[4, 3] <- 5
  spmat[1, 4] <- 6
  spmat[3, 4] <- 7

  list_ind_row <- list_ind_col <- list(2, 1:5, sample(5, replace = TRUE))

  for (compact in c(TRUE, FALSE)) {
    X <- as_SFBM(spmat, compact = compact)
    for (ind_row in list_ind_row) {
      for (ind_col in list_ind_col) {
        expect_identical(X$dense_acc(ind_row, ind_col),
                         as.matrix(spmat[ind_row, ind_col, drop = FALSE]))
      }
    }
  }

  spmat2 <- Matrix::rsparsematrix(10, 10, density = 0.5, symmetric = TRUE)
  diag(spmat2) <- 10
  X2 <- as_SFBM_corr_compact(Matrix::cov2cor(spmat2))
  for (ind_row in list_ind_row) {
    for (ind_col in list_ind_col) {
      expect_identical(X2$dense_acc(ind_row, ind_col),
                       as.matrix(X2[ind_row, ind_col]))
    }
  }
})

################################################################################
