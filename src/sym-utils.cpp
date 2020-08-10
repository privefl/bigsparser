/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

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
