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
#'
#' @importFrom bigassertr assert_class assert_noexist assert_dir
#' @importFrom methods new
#'
#' @exportClass SFBM
#'
SFBM_RC <- methods::setRefClass(

  "SFBM",

  fields = list(
    extptr = "externalptr",
    nrow = "numeric",
    ncol = "numeric",
    nval = "numeric",
    p = "integer",
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

    sbk = function() .self$backingfile,
    rds = function() sub("\\.sbk$", ".rds", .self$sbk),
    is_saved = function() file.exists(.self$rds)
  ),

  methods = list(
    initialize = function(spmat, backingfile) {

      sbkfile <- path.expand(paste0(backingfile, ".sbk"))

      assert_noexist(sbkfile)
      assert_dir(dirname(sbkfile))

      write_indval(sbkfile, spmat@i, spmat@x)

      .self$backingfile  <- normalizePath(sbkfile)
      .self$nrow         <- spmat@Dim[1]
      .self$ncol         <- spmat@Dim[2]
      .self$nval         <- length(spmat@x)
      .self$p            <- spmat@p
      .self$extptr       <- NIL_PTR

      .self
    },

    save = function() {
      saveRDS(.self, .self$rds)
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
SFBM_RC$lock("nrow")

################################################################################

#' Convert to SFBM
#'
#' Convert a dgCMatrix to an SFBM.
#'
#' @param spmat A dgCMatrix (non-symmetric sparse matrix of type 'double').
#' @param backingfile Path to file where to store data. Extension `.sbk` is
#'   automatically added.
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

  assert_class(spmat, "dgCMatrix")

  new("SFBM", spmat = spmat, backingfile = backingfile)
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
