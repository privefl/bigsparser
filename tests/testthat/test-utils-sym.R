################################################################################

test_that("col_count_sym() works", {

  replicate(100, {

    N <- 300
    spmat <- Matrix::rsparsematrix(N, N, 0.1, symmetric = TRUE)
    expect_identical(col_count_sym(spmat@p, spmat@i),
                     Matrix::colSums(spmat != 0))
  })

})

################################################################################
