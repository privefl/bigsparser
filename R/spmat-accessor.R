################################################################################

# Transform negative or boolean indices to positive indices
transform_ind <- function(k, lim) {

  if (missing(k))
    return(seq_len(lim))

  if (is.character(k))
    stop2("Character subsetting is not allowed.")

  res <- seq_len(lim)[k]

  if (anyNA(res))
    stop2("Error when subsetting (missing values? out of bounds?)")

  res
}


spmat_accessor <- function(compact, corr = FALSE) {

  function(x, i, j, ..., drop = FALSE) {

    if (!missing(drop))
      stop2("Parameter 'drop' is not implemented; do not use it.")

    nargs <- nargs()

    if (nargs == 2) {
      if (missing(i)) {
        nargs <- 3  # x[] is the same as x[,]
      } else {
        stop2("Vector subsetting is not allowed.")
      }
    }

    if (nargs == 3) {

      ind_row <- transform_ind(i, nrow(x))
      ind_col <- transform_ind(j, ncol(x))

      res <- if (corr) {
        access_subset_corr_compact(x, ind_row, ind_col)
      } else if (compact) {
        access_subset_compact(x, ind_row, ind_col)
      } else {
        access_subset(x, ind_row, ind_col)
      }

      sub <- new("dgCMatrix")
      sub@Dim <- c(length(ind_row), length(ind_col))
      sub@i <- res$i
      sub@p <- as.integer(res$p)
      sub@x <- res$x

      return(sub)
    }
  }
}

################################################################################

#' Accessor methods for class `SFBM`.
#'
#' @param x A [SFBM][SFBM-class] object.
#' @param i A vector of indices (or nothing). You can use positive and negative
#'   indices, and also logical indices (that are recycled).
#' @param j A vector of indices (or nothing). You can use positive and negative
#'   indices, and also logical indices (that are recycled).
#' @param ... Not used. Just to make [nargs] work.
#' @param drop Not implemented; always return a sparse matrix (`drop = FALSE`).
#'
#' @rdname crochet
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @export
#'
#' @examples
#' spmat <- Matrix::Diagonal(4, 0:3)
#' spmat[4, 2] <- 5
#' spmat[1, 4] <- 6
#' spmat[3, 4] <- 7
#' spmat
#'
#' X <- as_SFBM(spmat)
#' X[1:3, 2:3]
#' X[, 4]   # parameter drop is not implemented
#' X[-1, 3:4]
#'
#' X2 <- as_SFBM(spmat, compact = TRUE)
#' X2[1:3, 2:3]
#'
setMethod('[', signature(x = "SFBM"), spmat_accessor(compact = FALSE))


#' @rdname crochet
#'
#' @export
#'
setMethod('[', signature(x = "SFBM_compact"), spmat_accessor(compact = TRUE))


#' @rdname crochet
#'
#' @export
#'
setMethod('[', signature(x = "SFBM_corr_compact"), spmat_accessor(corr = TRUE))

################################################################################
