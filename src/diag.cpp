/******************************************************************************/

#include <bigsparser/SFBM.h>
#include <bigsparser/SFBM-corr-compact.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector diag_sfbm(Environment X) {

  XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  const double * data = sfbm->i_x();
  int n = std::min(sfbm->nrow(), sfbm->ncol());

  NumericVector diag(n);  // init with 0s

  for (int j = 0; j < n; j++) {

    size_t lo = 2 * p[j];
    size_t up = 2 * p[j + 1];

    for (size_t k = lo; k < up; k += 2) {
      int i = data[k];
      if (i >= j) {
        if (i == j) diag[j] = data[k + 1];
        break;
      }
    }
  }

  return diag;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector diag_sfbm_compact(Environment X) {

  XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  const IntegerVector first_i = X["first_i"];
  const double * data = sfbm->i_x();
  int n = std::min(sfbm->nrow(), sfbm->ncol());

  NumericVector diag(n);  // init with 0s

  for (int j = 0; j < n; j++) {

    int min_i = first_i[j];

    if (min_i >= 0) {

      int shift = j - min_i;
      if (shift >= 0) {
        size_t k = p[j] + shift;
        if (k < p[j + 1]) diag[j] = data[k];
      }
    }
  }

  return diag;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector diag_sfbm_corr_compact(Environment X) {

  XPtr<SFBM_corr_compact> sfbm = X["address"];
  const NumericVector p = X["p"];
  const IntegerVector first_i = X["first_i"];
  const int16_t * data = sfbm->i_x();
  int n = std::min(sfbm->nrow(), sfbm->ncol());

  NumericVector diag(n);  // init with 0s

  for (int j = 0; j < n; j++) {

    int min_i = first_i[j];

    if (min_i >= 0) {

      int shift = j - min_i;
      if (shift >= 0) {
        size_t k = p[j] + shift;
        if (k < p[j + 1]) diag[j] = data[k] / 32767.0;
      }
    }
  }

  return diag;
}

/******************************************************************************/
