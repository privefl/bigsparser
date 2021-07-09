/******************************************************************************/

#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
IntegerVector col_count_compact(const std::vector<size_t>& p,
                                const IntegerVector& i) {

  int m = p.size() - 1;
  IntegerVector res(m);

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    if (up == lo) {
      res[j] = 0;
    } else {
      res[j] = i[up - 1] - i[lo] + 1;
    }
  }

  return res;
}

/******************************************************************************/

// [[Rcpp::export]]
IntegerVector write_compact_val(const char * filename,
                                const std::vector<size_t>& p,
                                const IntegerVector& i,
                                const NumericVector& x) {

  std::ofstream myFile(filename, std::ios::out | std::ios::app | std::ios::binary);

  int m = p.size() - 1;
  IntegerVector first_i(m);

  double ZERO = 0, x_k = 0;
  const char * ptr_char_0 = reinterpret_cast<const char*>(&ZERO);
  const char * ptr_char_x = reinterpret_cast<const char*>(&x_k);

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    if (up == lo) {
      first_i[j] = -1;
    } else {
      double cur_i = first_i[j] = i[lo];
      for (size_t k = lo; k < up; k++) {
        while (cur_i < i[k]) {
          myFile.write(ptr_char_0, 8);
          cur_i++;
        }
        x_k = x[k];
        myFile.write(ptr_char_x, 8);
        cur_i++;
      }
    }
  }

  myFile.close();

  return first_i;
}

/******************************************************************************/

// [[Rcpp::export]]
ListOf<IntegerVector> col_range_sym(const std::vector<size_t>& p,
                                    const IntegerVector& i) {

  int m = p.size() - 1;
  IntegerVector first_i(m, -1), last_i(m, -1);

  for (int j = 0; j < m; j++) {

    size_t lo = p[j];
    size_t up = p[j + 1];

    if (up > lo) {
      first_i[j] = i[lo];
      for (size_t k = lo; k < up; k++)
        if (last_i[i[k]] < j) last_i[i[k]] = j;
    }
  }

  return List::create(first_i, last_i);
}

/******************************************************************************/
