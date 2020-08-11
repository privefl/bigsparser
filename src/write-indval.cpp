/******************************************************************************/

#include <rmio/create-file.hpp>
#include <bigsparser/SFBM.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
void write_indval(std::string filename,
                  const IntegerVector& i,
                  const NumericVector& x,
                  size_t offset_p,
                  int offset_i = 0) {

  size_t K = i.size();

  if (offset_p == 0) {
    create_file(filename.c_str(), K, 16);
  } else {
    append_file(filename.c_str(), K, 16);
  }

  mio::mmap_sink rw_mmap;
  std::error_code error;
  rw_mmap.map(filename, error);
  if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

  double * data = reinterpret_cast<double*>(rw_mmap.data());

  size_t k2 = 2 * offset_p;
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
                               size_t offset_p,
                               int offset_i = 0) {

  IntegerVector count = col_count_sym(p, i);
  int m = count.size();

  std::vector<size_t> data_offset(m);
  size_t K = offset_p;
  for (int j = 0; j < m; j++) {
    K += count[j];
    data_offset[j] = 2 * K;
  }

  if (offset_p == 0) {
    create_file(filename.c_str(), K, 16);
  } else {
    append_file(filename.c_str(), K - offset_p, 16);
  }

  mio::mmap_sink rw_mmap;
  std::error_code error;
  rw_mmap.map(filename, error);
  if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

  double * data = reinterpret_cast<double*>(rw_mmap.data());

  for (int j = m - 1; j >= 0; j--) {

    size_t lo = p[j];
    size_t up = p[j + 1];
    if (up == 0) continue;  // cannot have -1 with size_t

    for (size_t k = up - 1; k >= lo; k--) {

      int    ind = i[k];
      double val = x[k];

      // write (i, j, x)
      size_t where = data_offset[j];
      data[--where] = val;
      data[--where] = ind + offset_i;
      data_offset[j] = where;

      if (ind != j) {
        // write (j, i, x)
        where = data_offset[ind];
        data[--where] = val;
        data[--where] = j + offset_i;
        data_offset[ind] = where;
      }

      if (k == 0) break;  // cannot have -1 with size_t
    }
  }

  NumericVector new_p(m + 1);
  for (int j = 0; j < m; j++) new_p[j] = data_offset[j] / 2;
  new_p[m] = K;

  return new_p;
}

/******************************************************************************/
