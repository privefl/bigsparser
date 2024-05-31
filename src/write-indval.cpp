/******************************************************************************/

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
void write_indval(std::string filename,
                  const IntegerVector& i,
                  const NumericVector& x,
                  size_t offset_p,
                  int offset_i) {

  size_t K = x.size();

  mio::mmap_sink rw_mmap;
  std::error_code error;
  rw_mmap.map(filename, 16 * offset_p, 16 * K, error);
  if (error)
    Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

  double * data = reinterpret_cast<double*>(rw_mmap.data());

  size_t k2 = 0;
  for (size_t k = 0; k < K; k++) {
    data[k2++] = i[k] + offset_i;
    data[k2++] = x[k];
  }
}

/******************************************************************************/

// [[Rcpp::export]]
IntegerVector col_count_sym(std::vector<size_t> p,
                            const IntegerVector& i) {

  int m = p.size() - 1;
  IntegerVector res(m);

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    for (size_t k = lo; k < up; k++) {
      (res[j])++;
      if (i[k] != j) (res[i[k]])++;
    }
  }

  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector write_indval_sym(std::string filename,
                               std::vector<size_t> p,
                               const IntegerVector& i,
                               const NumericVector& x,
                               const IntegerVector& col_count,
                               size_t offset_p,
                               int offset_i) {

  int m = col_count.size();

  std::vector<size_t> data_offset(m);
  size_t K = 0;
  for (int j = 0; j < m; j++) {
    K += col_count[j];
    data_offset[j] = 2 * K;
  }

  mio::mmap_sink rw_mmap;
  std::error_code error;
  rw_mmap.map(filename, 16 * offset_p, 16 * K, error);
  if (error)
    Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

  double * data = reinterpret_cast<double*>(rw_mmap.data());

  for (int j = m - 1; j >= 0; j--) {

    size_t lo = p[j];
    size_t up = p[j + 1];
    if (up == 0) continue;  // cannot have -1 with size_t

    for (size_t k = up - 1; k >= lo; k--) {

      int    ind = i[k];
      double val = x[k];

      // write (i, j, x)
      size_t where1 = data_offset[j];
      data[--where1] = val;
      data[--where1] = ind + offset_i;
      data_offset[j] = where1;

      if (ind != j) {
        // write (j, i, x)
        size_t where2 = data_offset[ind];
        data[--where2] = val;
        data[--where2] = j + offset_i;
        data_offset[ind] = where2;
      }

      if (k == 0) break;  // cannot have -1 with size_t
    }
  }

  NumericVector new_p(m + 1);
  size_t K2 = 0;
  new_p[0] = offset_p;
  for (int j = 0; j < m; j++) {
    // data_offset has been reduced such that (data_offset[j + 1] -> data_offset[j])
    if (data_offset[j] != (2 * K2)) Rcpp::stop("This is a bug.");
    K2 += col_count[j];
    new_p[j + 1] = new_p[j] + col_count[j];
  }

  return new_p;
}

/******************************************************************************/
