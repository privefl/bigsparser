################################################################################

NIL_PTR <- methods::new("externalptr")

assert_sparse_matrix <- function(spmat) {

  if (inherits(spmat, "dgCMatrix")) {
    FALSE
  } else if (inherits(spmat, "dsCMatrix")) {
    if (!identical(spmat@uplo, "U"))
      stop2("Only upper-triangle symmetric matrices are implemented.")
    TRUE
  } else {
    stop2("Only classes 'dgCMatrix' and 'dsCMatrix' are currently implemented.")
  }
}

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
    initialize = function(spmat, backingfile) {

      symmetric <- assert_sparse_matrix(spmat)

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

      symmetric <- assert_sparse_matrix(spmat)

      assert_one_int(offset_i)
      assert_pos(offset_i, strict = FALSE)

      sbkfile <- .self$sbk
      assert_exist(sbkfile)

      offset_p <- tail(.self$p, 1)

      if (symmetric) {
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
#' @param compact Whether to use a compact format? Default is `FALSE`.
#'   This is useful when non-zero values in columns are contiguous (or almost).
#'
#' @return The new [SFBM][SFBM-class].
#'
#' @rdname SFBM-class
#' @export
#'
#' @example examples/example-SFBM.R
#'
as_SFBM <- function(spmat, backingfile = tempfile(), compact = FALSE) {

  if (compact) {
    new("SFBM_compact", spmat = spmat, backingfile = backingfile)
  } else {
    new("SFBM",         spmat = spmat, backingfile = backingfile)
  }
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

#' Class SFBM_compact
#'
#' A reference class for storing and accessing sparse matrix-like data stored
#' in files on disk, in a compact format (when non-zero values in columns are
#' contiguous).
#'
#' @details
#' It inherits the fields and methods from class [SFBM][SFBM-class].
#'
#' @importFrom bigassertr assert_exist assert_noexist assert_dir
#' @importFrom bigassertr assert_class assert_pos assert_one_int stop2
#' @importFrom methods new
#' @importFrom utils head tail
#'
#' @exportClass SFBM_compact
#'
SFBM_compact_RC <- methods::setRefClass(

  "SFBM_compact",

  contains = "SFBM",

  fields = list(
    first_i = "integer",

    #### Active bindings
    address = function() {

      if (identical(.self$extptr, NIL_PTR)) {
        .self$extptr <- getXPtrSFBM_compact(
          path    = .self$backingfile,
          n       = .self$nrow,
          m       = .self$ncol,
          p       = .self$p,
          first_i = .self$first_i
        )
      }

      .self$extptr
    }
  ),

  methods = list(
    initialize = function(spmat, backingfile) {

      symmetric <- assert_sparse_matrix(spmat)

      sbkfile <- path.expand(paste0(backingfile, ".sbk"))
      assert_noexist(sbkfile)
      assert_dir(dirname(sbkfile))

      res <- write_val_compact(sbkfile, spmat@p, spmat@i, spmat@x, 0, 0L,
                               symmetric = symmetric)

      .self$backingfile <- normalizePath(sbkfile)
      .self$nrow        <- spmat@Dim[1]
      .self$first_i     <- res[[1]]
      .self$p           <- res[[2]]
      .self$extptr      <- NIL_PTR

      .self
    },

    add_columns = function(spmat, offset_i) {

      symmetric <- assert_sparse_matrix(spmat)

      assert_one_int(offset_i)
      assert_pos(offset_i, strict = FALSE)

      sbkfile <- .self$sbk
      assert_exist(sbkfile)

      offset_p <- tail(.self$p, 1)

      res <- write_val_compact(sbkfile, spmat@p, spmat@i, spmat@x,
                               offset_p, offset_i, symmetric = symmetric)

      .self$nrow    <- max(.self$nrow, spmat@Dim[1] + offset_i)
      .self$first_i <- c(.self$first_i, res[[1]])
      .self$p       <- c(head(.self$p, -1), res[[2]])
      .self$extptr  <- NIL_PTR

      .self
    },

    show = function() {  # TODO: add type when implementing 'float'
      cat(sprintf(
        "A compact Sparse Filebacked Big Matrix with %s rows and %s columns.\n",
        .self$nrow, .self$ncol))
      invisible(.self)
    }
  )
)

################################################################################
