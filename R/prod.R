################################################################################

#' Products with a vector
#'
#' Products between an [SFBM][SFBM-class] and a vector.
#'
#' @param X An [SFBM][SFBM-class].
#' @param y A vector of same size of the number of columns of `X` for
#'   `sp_prodVec()` abd as the number of rows of `X` for `sp_cprodVec()`.
#'
#' @return
#' - `sp_prodVec()`: the vector which is equivalent to `X %*% y`
#'   if `X` was a dgCMatrix.
#' - `sp_cprodVec()`: the vector which is equivalent to `Matrix::crossprod(X, y)`
#'   if `X` was a dgCMatrix.
#'
#' @export
#'
#' @importFrom bigassertr assert_lengths
#'
#' @examples
#' spmat <- Matrix::rsparsematrix(1000, 1000, 0.01)
#' X <- as_SFBM(spmat)
#' sp_prodVec(X, rep(1, 1000))
#' sp_cprodVec(X, rep(1, 1000))
#'
sp_prodVec <- function(X, y) {
  assert_lengths(seq_len(ncol(X)), y)
  prodVec(X, y)
}

################################################################################

#' @rdname sp_prodVec
#' @export
sp_cprodVec <- function(X, y) {
  assert_lengths(seq_len(nrow(X)), y)
  cprodVec(X, y)
}

################################################################################
