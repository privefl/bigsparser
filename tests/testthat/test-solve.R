################################################################################

test_that("sp_solve_sym() works", {

  skip_on_cran()

  replicate(20, {

    N <- 300
    spmat <- Matrix::rsparsematrix(N, N, 0.01, symmetric = TRUE) + Matrix::Diagonal(N, 1e-2)
    X <- bigsparser::as_SFBM(as(spmat, "dgCMatrix"))
    b <- runif(N)

    test1 <- try(as.vector(Matrix::solve(spmat, b)))
    if (class(test1) != "try-error") {
      test2 <- sp_solve_sym(X, b)
      expect_equal(test1, test2)
    }

    test3 <- as.vector(Matrix::solve(spmat + Matrix::Diagonal(N, 1:N), b))
    test4 <- sp_solve_sym(X, b, add_to_diag = 1:N)
    expect_equal(test3, test4)

    test5 <- as.vector(Matrix::solve(spmat + Matrix::Diagonal(N), b))
    test6 <- sp_solve_sym(X, b, add_to_diag = 1)
    expect_equal(test5, test6)
  })

})

################################################################################

test_that("sp_solve_sym() works with compact SFBM", {

  skip_on_cran()

  replicate(20, {

    N <- 300
    spmat <- Matrix::rsparsematrix(N, N, 0.01, symmetric = TRUE) + Matrix::Diagonal(N, 1e-2)
    X <- bigsparser::as_SFBM(spmat, compact = TRUE)
    b <- runif(N)

    test1 <- try(as.vector(Matrix::solve(spmat, b)))
    if (class(test1) != "try-error") {
      test2 <- sp_solve_sym(X, b)
      expect_equal(test1, test2)
    }

    test3 <- as.vector(Matrix::solve(spmat + Matrix::Diagonal(N, 1:N), b))
    test4 <- sp_solve_sym(X, b, add_to_diag = 1:N)
    expect_equal(test3, test4)

    test5 <- as.vector(Matrix::solve(spmat + Matrix::Diagonal(N), b))
    test6 <- sp_solve_sym(X, b, add_to_diag = 1)
    expect_equal(test5, test6)
  })

})

################################################################################
