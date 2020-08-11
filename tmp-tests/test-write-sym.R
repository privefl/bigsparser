library(Matrix)
spmat <- rsparsematrix(5, 5, nnz = 10, symmetric = TRUE)

Rcpp::sourceCpp('src/write-indval.cpp')
col_count_sym(spmat@p, spmat@i)
Matrix::colSums(spmat != 0)

tmp <- tempfile()
test <- write_indval_sym(tmp, spmat@p, spmat@i, spmat@x, 0)


spmat@p
spmat2 <- as(spmat, "dgCMatrix"); spmat2@p

true <- matrix(readBin(bigsparser::as_SFBM(spmat2)$sbk,
                       what = 1, n = 100), nrow = 2)
rbind(matrix(readBin(tmp, what = 1, n = 100), nrow = 2), true)

all.equal(readBin(bigsparser::as_SFBM(spmat2)$sbk, what = 1, n = 100),
          readBin(tmp, what = 1, n = 100))


write_indval_sym_r <- function(spmat, offset = 0) {

  p <- spmat@p
  i <- spmat@i
  x <- spmat@x

  count = col_count_sym(p, i);
  m = length(count);

  data_offset <- c()
  K = offset;
  for (j in 1:m) {
    K = K +  count[j];
    data_offset[j] = 2 * K;
  }
  c(K, sum(count))

  data <- bigstatsr::FBM(2, K)

  for (j in m:1) {

    # j <- m

    lo = p[j];
    up = p[j + 1];

    if (up > lo) for (k in seq(up, lo + 1)) {

      # k <- up

      ind = i[k] + 1L;
      val = x[k];

      # write (i, j, x)
      where = data_offset[j];
      print(where)
      data[where] = val;
      where <- where - 1L
      data[where] = ind - 1L;
      where <- where - 1L
      data_offset[j] = where;

      print(data_offset[])
      print(rbind(data[], true))

      # write (j, i, x)
      if (ind != j) {
        where = data_offset[ind];
        print(where)
        data[where] = val;
        where <- where - 1L
        data[where] = j - 1L;
        where <- where - 1L
        data_offset[ind] = where;

        print(data_offset[])
        print(rbind(data[], true))
      }
    }
  }
}

# debug(write_indval_sym_r)
write_indval_sym_r(spmat)

true <- bigsparser::as_SFBM(as(spmat, "dgCMatrix"), backingfile = tmp)
readBin(true$sbk, what = 1, n = 1000)
all.equal(readBin(tmp,      what = 1, n = 1000),
          readBin(true$sbk, what = 1, n = 1000))


sfbm <- as_SFBM(spmat)
dim(sfbm)
sfbm$add_columns(spmat, 0)
dim(sfbm)
sfbm$add_columns(spmat, 10)
