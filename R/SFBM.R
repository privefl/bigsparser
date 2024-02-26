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
#'   - `$add_columns()`: Add new columns from a 'dgCMatrix' or a 'dsCMatrix'.
#'   - `$dense_acc()`: Equivalent to `as.matrix(.[ind_row, ind_col])`. Use with
#'     caution; `ind_row` and `ind_col` must be positive indices within range.
#'
#' @importFrom bigassertr assert_class assert_pos assert_one_int stop2
#' @importFrom methods new
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
    nval = function() utils::tail(.self$p, 1),

    sbk = function() .self$backingfile,
    rds = function() sub("\\.sbk$", ".rds", .self$sbk),
    is_saved = function() file.exists(.self$rds)
  ),

  methods = list(
    initialize = function(spmat, backingfile) {

      symmetric <- assert_sparse_matrix(spmat)
      sbkfile <- paste0(backingfile, ".sbk")

      if (symmetric) {
        col_count <- col_count_sym(spmat@p, spmat@i)
        sbkfile <- rmio::file_create(sbkfile, 16 * sum(col_count))
        .self$p <- write_indval_sym(sbkfile, spmat@p, spmat@i, spmat@x,
                                    col_count, offset_p = 0L, offset_i = 0L)
      } else {
        sbkfile <- rmio::file_create(sbkfile, 16 * length(spmat@x))
        write_indval(sbkfile, spmat@i, spmat@x, offset_p = 0L, offset_i = 0L)
        .self$p <- spmat@p
      }

      .self$backingfile <- normalizePath(sbkfile)
      .self$nrow        <- spmat@Dim[1]
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

      offset_p <- .self$nval

      ## reset pointers -> need this before resizing
      .self$extptr <- NIL_PTR
      gc()

      if (symmetric) {
        col_count <- col_count_sym(spmat@p, spmat@i)
        sbkfile <- rmio::file_resize_off(.self$sbk, 16 * sum(col_count))
        new_p <- write_indval_sym(sbkfile, spmat@p, spmat@i, spmat@x,
                                  col_count, offset_p, offset_i)
      } else {
        sbkfile <- rmio::file_resize_off(.self$sbk, 16 * length(spmat@x))
        write_indval(sbkfile, spmat@i, spmat@x, offset_p, offset_i)
        new_p <- spmat@p + as.double(offset_p)
      }

      .self$nrow <- max(.self$nrow, spmat@Dim[1] + offset_i)
      .self$p    <- c(.self$p, new_p[-1])

      .self
    },

    dense_acc = function(ind_row, ind_col) {
      access_dense_subset(.self, ind_row, ind_col)
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
#' Convert a 'dgCMatrix' or 'dsCMatrix' to an SFBM.
#'
#' @param spmat A 'dgCMatrix' (non-symmetric sparse matrix of type 'double')
#'   or 'dsCMatrix' (symmetric sparse matrix of type 'double').
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
    methods::new("SFBM_compact", spmat = spmat, backingfile = backingfile)
  } else {
    methods::new("SFBM",         spmat = spmat, backingfile = backingfile)
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

      col_range <- `if`(symmetric, range_col_sym, range_col)(spmat@p, spmat@i)
      first_i_  <- col_range[[1]]
      col_count <- col_range[[2]] - first_i_ + 1L

      sbkfile <- rmio::file_create(paste0(backingfile, ".sbk"), 8 * sum(col_count))

      .self$p <- write_val_compact(sbkfile, spmat@p, spmat@i, spmat@x,
                                   first_i_, col_count, offset_p = 0L, symmetric)

      .self$backingfile <- normalizePath(sbkfile)
      .self$nrow        <- spmat@Dim[1]
      .self$first_i     <- first_i_
      .self$extptr      <- NIL_PTR

      .self
    },

    add_columns = function(spmat, offset_i) {

      symmetric <- assert_sparse_matrix(spmat)

      assert_one_int(offset_i)
      assert_pos(offset_i, strict = FALSE)

      offset_p <- .self$nval

      col_range <- `if`(symmetric, range_col_sym, range_col)(spmat@p, spmat@i)
      first_i_  <- col_range[[1]]
      col_count <- col_range[[2]] - first_i_ + 1L

      ## reset pointers -> need this before resizing
      .self$extptr <- NIL_PTR
      gc()

      sbkfile <- rmio::file_resize_off(.self$sbk, 8 * sum(col_count))

      new_p <- write_val_compact(sbkfile, spmat@p, spmat@i, spmat@x,
                                first_i_, col_count, offset_p, symmetric)

      .self$nrow    <- max(.self$nrow, spmat@Dim[1] + offset_i)
      .self$p       <- c(.self$p, new_p[-1])
      .self$first_i <- c(.self$first_i,
                         ifelse(first_i_ >= 0, first_i_ + as.integer(offset_i), first_i_))

      .self
    },

    dense_acc = function(ind_row, ind_col) {
      access_dense_subset_compact(.self, ind_row, ind_col)
    },

    show = function() {
      cat(sprintf(
        "A compact Sparse Filebacked Big Matrix with %s rows and %s columns.\n",
        .self$nrow, .self$ncol))
      invisible(.self)
    }
  )
)

################################################################################

#' Class SFBM_corr_compact
#'
#' A reference class for storing and accessing from disk a sparse correlation
#' matrix where non-zero values in columns are mostly contiguous. It rounds
#' correlation values with precision 1/32767 to store them using 2 bytes only.
#' This class has been specifically designed for package 'bigsnpr'.
#'
#' @details
#' It inherits the fields and methods from class [SFBM_compact][SFBM-class].
#'
#' @exportClass SFBM_corr_compact
#'
SFBM_corr_compact_RC <- methods::setRefClass(

  "SFBM_corr_compact",

  contains = "SFBM_compact",

  fields = list(

    #### Active bindings
    address = function() {

      if (identical(.self$extptr, NIL_PTR)) {
        .self$extptr <- getXPtrSFBM_corr_compact(
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

      col_range <- `if`(symmetric, range_col_sym, range_col)(spmat@p, spmat@i)
      first_i_  <- col_range[[1]]
      col_count <- col_range[[2]] - first_i_ + 1L

      sbkfile <- rmio::file_create(paste0(backingfile, ".sbk"), 2 * sum(col_count))

      .self$p <- write_val_corr_compact(
        sbkfile, spmat@p, spmat@i, spmat@x,
        first_i_, col_count, offset_p = 0L, symmetric)

      .self$backingfile <- normalizePath(sbkfile)
      .self$nrow        <- spmat@Dim[1]
      .self$first_i     <- first_i_
      .self$extptr      <- NIL_PTR

      .self
    },

    add_columns = function(spmat, offset_i) {

      symmetric <- assert_sparse_matrix(spmat)

      assert_one_int(offset_i)
      assert_pos(offset_i, strict = FALSE)

      offset_p <- .self$nval

      col_range <- `if`(symmetric, range_col_sym, range_col)(spmat@p, spmat@i)
      first_i_  <- col_range[[1]]
      col_count <- col_range[[2]] - first_i_ + 1L

      ## reset pointers -> need this before resizing
      .self$extptr <- NIL_PTR
      gc()

      sbkfile <- rmio::file_resize_off(.self$sbk, 2 * sum(col_count))

      new_p <- write_val_corr_compact(
        sbkfile, spmat@p, spmat@i, spmat@x,
        first_i_, col_count, offset_p, symmetric)

      .self$nrow    <- max(.self$nrow, spmat@Dim[1] + offset_i)
      .self$p       <- c(.self$p, new_p[-1])
      .self$first_i <- c(.self$first_i,
                         ifelse(first_i_ >= 0, first_i_ + as.integer(offset_i), first_i_))

      .self
    },

    dense_acc = function(ind_row, ind_col) {
      access_dense_subset_corr_compact(.self, ind_row, ind_col)
    },

    show = function() {
      cat(sprintf(
        "A compact SFBM correlation with %s rows and %s columns.\n",
        .self$nrow, .self$ncol))
      invisible(.self)
    }
  )
)

################################################################################

#' Convert to SFBM_corr_compact
#'
#' Convert a 'dgCMatrix' or 'dsCMatrix' to an SFBM_corr_compact.
#'
#' @param spmat A 'dgCMatrix' (non-symmetric sparse matrix of type 'double')
#'   or 'dsCMatrix' (symmetric sparse matrix of type 'double').
#' @param backingfile Path to file where to store data. Extension `.sbk` is
#'   automatically added.
#'
#' @return The new [SFBM_corr_compact][SFBM_corr_compact-class].
#'
#' @rdname SFBM_corr_compact-class
#' @export
#'
#' @examples
#' spmat2 <- as(cor(iris[1:4]), "dsCMatrix")
#' (X2 <- as_SFBM_corr_compact(spmat2))
#' (bin <- readBin(X2$sbk, what = integer(), size = 2, n = 100))
#' matrix(bin / 32767, 4)
#' spmat2
#'
as_SFBM_corr_compact <- function(spmat, backingfile = tempfile()) {
  methods::new("SFBM_corr_compact", spmat = spmat, backingfile = backingfile)
}

################################################################################
