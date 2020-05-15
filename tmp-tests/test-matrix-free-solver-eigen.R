N <- 100
spmat <- Matrix::rsparsematrix(N, N, 0.01, symmetric = TRUE) + Matrix::Diagonal(N, 1)
X <- bigsparser::as_SFBM(as(spmat, "dgCMatrix"))
b <- runif(N)
system.time(
  test <- Matrix::solve(spmat, b)
) # 34 sec

Rcpp::sourceCpp('tmp-tests/test-matrix-free-solver-eigen3.cpp')
system.time(
  test2 <- spsolve(X, b, use_CG = FALSE)
) # 5 sec for CG + MINRES

plot(test, test2)
