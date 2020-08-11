################################################################################

NIL_PTR <- methods::new("externalptr")

################################################################################

#' Class SFBM
#'
#' A reference class for storing and accessing sparse matrix-like data stored
#' in files on disk.
#'
#' @details
#' An object of class SFBM has many fields:
#'   - `$address`: address of the external pointer containing the underlying
#'     C++ object to be used as a `XPtr<SFBM>` in C++ code
#'   - `$extptr`: (internal) use `$address` instead
#'   - `$nrow`: number of rows
#'   - `$ncol`: number of columns
#'   - `$nval`: number of non-zero values
#'   - `$p`: vector of column positions
#'   - `$backingfile` or `$sbk`: File with extension 'sbk' that stores the
#'     data of the SFBM
#'   - `$rds`: 'rds' file (that may not exist) corresponding to the 'sbk' file
#'   - `$is_saved`: whether this object is stored in `$rds`?
#'
#' And some methods:
#'   - `$save()`: Save the SFBM object in `$rds`. Returns the SFBM.
#'   - `$add_columns()`: Add new columns from another sparse dgCMatrix.
#'
#' @importFrom bigassertr assert_exist assert_noexist assert_dir
#' @importFrom bigassertr assert_class assert_pos assert_one_int stop2
#' @importFrom methods new
#' @importFrom utils head tail
#'
#' @exportClass SFBM
#'
SFBM_RC <- methods::setRefClass(

  "SFBM",

  fields = list(
    extptr = "externalptr",
    nrow = "numeric",
    p = "numeric",
    backingfile = "character",

    #### Active bindings
    address = function() {

      if (identical(.self$extptr, NIL_PTR)) {
        .self$extptr <- getXPtrSFBM(
          path = .self$backingfile,
          n    = .self$nrow,
          m    = .self$ncol,
          p    = .self$p
        )
      }

      .self$extptr
    },

    ncol = function() length(.self$p) - 1L,
    nval = function() tail(.self$p, 1),

    sbk = function() .self$backingfile,
    rds = function() sub("\\.sbk$", ".rds", .self$sbk),
    is_saved = function() file.exists(.self$rds)
  ),

  methods = list(
    initialize = function(spmat, backingfile, symmetric = FALSE) {

      sbkfile <- path.expand(paste0(backingfile, ".sbk"))
      assert_noexist(sbkfile)
      assert_dir(dirname(sbkfile))

      if (symmetric) {
        new_p <- write_indval_sym(sbkfile, spmat@p, spmat@i, spmat@x, 0L)
      } else {
        write_indval(sbkfile, spmat@i, spmat@x, 0L)
      }

      .self$backingfile <- normalizePath(sbkfile)
      .self$nrow        <- spmat@Dim[1]
      .self$p           <- `if`(symmetric, new_p, spmat@p)
      .self$extptr      <- NIL_PTR

      .self
    },

    save = function() {
      saveRDS(.self, .self$rds)
      .self
    },

    add_columns = function(spmat, offset_i) {

      if (!inherits(spmat, "dgCMatrix") && !inherits(spmat, "dsCMatrix"))
        stop2("Only classes 'dgCMatrix' and 'dsCMatrix' are currently implemented.")

      assert_one_int(offset_i)
      assert_pos(offset_i, strict = FALSE)

      sbkfile <- .self$sbk
      assert_exist(sbkfile)

      offset_p <- tail(.self$p, 1)

      if (inherits(spmat, "dsCMatrix")) {
        new_p <- write_indval_sym(sbkfile, spmat@p, spmat@i, spmat@x,
                                  offset_p, offset_i)
      } else {
        write_indval(sbkfile, spmat@i, spmat@x, offset_p, offset_i)
        new_p <- spmat@p + 0 + offset_p
      }

      .self$nrow   <- max(.self$nrow, spmat@Dim[1] + offset_i)
      .self$p      <- c(head(.self$p, -1), new_p)
      .self$extptr <- NIL_PTR

      .self
    },

    show = function() {
      cat(sprintf(
        "A Sparse Filebacked Big Matrix with %s rows and %s columns.\n",
        .self$nrow, .self$ncol))
      invisible(.self)
    }
  )
)

################################################################################

#' Convert to SFBM
#'
#' Convert a dgCMatrix or dsCMatrix to an SFBM.
#'
#' @param spmat A dgCMatrix (non-symmetric sparse matrix of type 'double')
#'   or dsCMatrix (symmetric sparse matrix of type 'double').
#' @param backingfile Path to file where to store data. Extension `.sbk` is
#'   automatically added.
#'
#' @return The new [SFBM][SFBM-class].
#'
#' @rdname SFBM-class
#' @export
#'
#' @examples
#' spmat <- Matrix::rsparsematrix(1000, 1000, 0.01)
#' class(spmat)
#' (X <- as_SFBM(spmat))
#'
as_SFBM <- function(spmat, backingfile = tempfile()) {

  if (inherits(spmat, "dgCMatrix")) {
    sym <- FALSE
  } else if (inherits(spmat, "dsCMatrix")) {
    sym <- TRUE
    if (!identical(spmat@uplo, "U"))
      stop2("Only upper-triangle symmetric matrices are implemented.")
  } else {
    stop2("Only classes 'dgCMatrix' and 'dsCMatrix' are currently implemented.")
  }

  new("SFBM", spmat = spmat, backingfile = backingfile, symmetric = sym)
}

################################################################################

#' Dimension and type methods for class `SFBM`.
#'
#' @param x An object of class [SFBM][SFBM-class].
#'
#' @rdname SFBM-methods
#' @export
setMethod("dim",    signature(x = "SFBM"), function(x) c(x$nrow, x$ncol))

#' @rdname SFBM-methods
#' @export
setMethod("length", signature(x = "SFBM"), function(x) prod(dim(x)))

################################################################################
