/******************************************************************************/

#include <rmio/create-file.hpp>
#include <bigsparser/SFBM.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

/******************************************************************************/

size_t sum_int_no_overflow(const IntegerVector& x) {

  size_t sum = 0;
  for (const int& x_k : x) sum += x_k;

  return sum;
}

/******************************************************************************/

// [[Rcpp::export]]
ListOf<IntegerVector> range_col(const std::vector<size_t>& p,
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
ListOf<IntegerVector> range_col_sym(const std::vector<size_t>& p,
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
        if (last_i[i[k]] < j) last_i[i[k]] = j;
      }
    }
  }

  return List::create(first_i, last_i);
}

/******************************************************************************/

// [[Rcpp::export]]
List write_val_compact(std::string filename,
                       const std::vector<size_t>& p,
                       const IntegerVector& i,
                       const NumericVector& x,
                       size_t offset_p,
                       int offset_i,
                       bool symmetric) {

  ListOf<IntegerVector> col_range =
    symmetric ? range_col_sym(p, i) : range_col(p, i);

  IntegerVector first_i = col_range[0];
  IntegerVector col_count = col_range[1] - first_i + 1.0;
  if (is_true(any(col_count < 0))) Rcpp::stop("This is a bug.");

  size_t K = sum_int_no_overflow(col_count);

  if (offset_p == 0) {
    create_file(filename.c_str(), K, 8);
  } else {
    append_file(filename.c_str(), K, 8);
  }

  mio::mmap_sink rw_mmap;
  std::error_code error;
  rw_mmap.map(filename, error);
  if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

  double * data = reinterpret_cast<double*>(rw_mmap.data());

  int m = first_i.size();
  std::vector<size_t> new_p(m + 1);
  new_p[0] = offset_p;
  for (int j = 0; j < m; j++) new_p[j + 1] = new_p[j] + col_count[j];

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    for (size_t k = lo; k < up; k++) {
      // write x for (i, j)
      size_t where = new_p[j] + (i[k] - first_i[j]);
      data[where] = x[k];
      if (symmetric) {
        // write x for (j, i)
        where = new_p[i[k]] + (j - first_i[i[k]]);
        data[where] = x[k];
      }
    }
  }

  for (int j = 0; j < m; j++)
    if (first_i[j] >= 0) first_i[j] += offset_i;

  return List::create(first_i, new_p);
}

/******************************************************************************/
