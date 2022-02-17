/******************************************************************************/

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
ListOf<IntegerVector> range_col(std::vector<size_t> p,
                                const IntegerVector& i) {

  int m = p.size() - 1;
  IntegerVector first_i(m, -1), last_i(m, -2);

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    if (up > lo) {
      first_i[j] = i[lo];
      last_i [j] = i[up - 1];
    }
  }

  return List::create(first_i, last_i);
}


// [[Rcpp::export]]
ListOf<IntegerVector> range_col_sym(std::vector<size_t> p,
                                    const IntegerVector& i) {

  int m = p.size() - 1;
  IntegerVector first_i(m, -1), last_i(m, -2);

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    if (up > lo) {
      first_i[j] = i[lo];
      double last = i[up - 1];
      if (last_i[j] < last) last_i[j] = last;
      for (size_t k = lo; k < up; k++) {
        if (first_i[i[k]] < 0) first_i[i[k]] = j;
        if ( last_i[i[k]] < j)  last_i[i[k]] = j;
      }
    }
  }

  return List::create(first_i, last_i);
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector write_val_compact(std::string filename,
                                std::vector<size_t> p,
                                const IntegerVector& i,
                                const NumericVector& x,
                                const IntegerVector& first_i,
                                const IntegerVector& col_count,
                                size_t offset_p,
                                bool symmetric) {

  if (is_true(any(col_count < 0)))
    Rcpp::stop("This is a bug.");

  int m = col_count.size();

  std::vector<size_t> data_offset(m);
  size_t K = 0;
  for (int j = 0; j < m; j++) {
    data_offset[j] = K;
    K += col_count[j];
  }

  mio::mmap_sink rw_mmap;
  std::error_code error;
  rw_mmap.map(filename, 8 * offset_p, 8 * K, error);
  if (error)
    Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

  double * data = reinterpret_cast<double*>(rw_mmap.data());

  // make sure holes are filled up with 0s
  for (size_t k = 0; k < K; k++) data[k] = 0;

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    for (size_t k = lo; k < up; k++) {
      // write x for (i, j)
      size_t where1 = data_offset[j] + (i[k] - first_i[j]);
      data[where1] = x[k];
      if (symmetric) {
        // write x for (j, i)
        size_t where2 = data_offset[i[k]] + (j - first_i[i[k]]);
        data[where2] = x[k];
      }
    }
  }

  NumericVector new_p(m + 1);
  new_p[0] = offset_p;
  for (int j = 0; j < m; j++) new_p[j + 1] = new_p[j] + col_count[j];

  return new_p;
}

/******************************************************************************/
