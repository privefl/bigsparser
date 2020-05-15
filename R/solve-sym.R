################################################################################

#' Solver for symmetric SFBM
#'
#' Solve Ax=b where A is a symmetric SFBM, and b is a vector.
#'
#' @param A A symmetric [SFBM][SFBM-class].
#' @param b A vector.
#' @param add_to_diag Vector (or single value) to *virtually* add to
#'   the diagonal of `A`. Default is 0s.
#' @param tol Tolerance for convergence. Default is `1e-10`.
#' @param maxiter Maximum number of iterations for convergence.
#'
#' @return The vector x, solution of Ax=b.
#' @export
#'
#' @examples
#' N <- 100
#' spmat <- Matrix::rsparsematrix(N, N, 0.01, symmetric = TRUE)
#' X <- bigsparser::as_SFBM(as(spmat, "dgCMatrix"))
#' b <- runif(N)
#'
#' test <- tryCatch(as.vector(Matrix::solve(spmat, b)), error = function(e) print(e))
#' test2 <- tryCatch(sp_solve_sym(X, b), error = function(e) print(e))
#'
#' test3 <- as.vector(Matrix::solve(spmat + Matrix::Diagonal(N, 1:N), b))
#' test4 <- sp_solve_sym(X, b, add_to_diag = 1:N)
#' all.equal(test3, test4)
#'
sp_solve_sym <- function(A, b,
                         add_to_diag = rep(0, ncol(A)),
                         tol = 1e-10,
                         maxiter = 10 * ncol(A)) {

  assert_class(A, "SFBM")

  if (length(add_to_diag) == 1)
    add_to_diag <- rep(add_to_diag, ncol(A))

  assert_lengths(seq_len(nrow(A)), seq_len(ncol(A)), add_to_diag, b)

  sp_solve_sym_eigen(A, b, add_to_diag, tol, maxiter)
}

################################################################################
