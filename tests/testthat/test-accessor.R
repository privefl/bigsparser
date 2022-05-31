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
    list(2, 1:5, c(TRUE, FALSE, TRUE, FALSE, TRUE), -3, -(3:5), rep(1, 5))

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
})

################################################################################
