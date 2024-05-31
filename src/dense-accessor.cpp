/******************************************************************************/

#include <bigsparser/SFBM.h>
#include <bigsparser/SFBM-corr-compact.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix access_dense_subset(Environment X,
                                  const IntegerVector& ind_row,
                                  const IntegerVector& ind_col) {

  XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  const double * data = sfbm->i_x();
  const IntegerVector ind_row2 = ind_row - 1;
  int N = sfbm->nrow();
  NumericVector cur_col(N);

  int n = ind_row.size();
  int m = ind_col.size();
  NumericMatrix mat(n, m);

  for (int j = 0; j < m; j++) {

    int j2 = ind_col[j] - 1;

    size_t lo = 2 * p[j2];
    size_t up = 2 * p[j2 + 1];

    // first pass: store all elements from column
    for (size_t k = lo; k < up; k += 2) {
      int    ind = data[k];
      double val = data[k + 1];
      cur_col[ind] = val;
    }

    // second pass: read requested indices, and store corresponding values
    for (int i = 0; i < n; i++)
      mat(i, j) = cur_col[ind_row2[i]];

    // third pass: reset column
    for (size_t k = lo; k < up; k += 2) {
      int ind = data[k];
      cur_col[ind] = 0;
    }
  }

  return mat;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix access_dense_subset_compact(Environment X,
                                          const IntegerVector& ind_row,
                                          const IntegerVector& ind_col) {

  XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  const IntegerVector first_i = X["first_i"];
  const double * data = sfbm->i_x();
  const IntegerVector ind_row2 = ind_row - 1;

  int n = ind_row.size();
  int m = ind_col.size();
  NumericMatrix mat(n, m);

  for (int j = 0; j < m; j++) {

    int j2 = ind_col[j] - 1;

    int min_i = first_i[j2];
    if (min_i >= 0) {

      size_t lo = p[j2];
      int nb = p[j2 + 1] - lo;

      for (int i = 0; i < n; i++) {

        int i2 = ind_row2[i];
        if (i2 < min_i) continue;  // before range
        int shift = i2 - min_i;
        if (shift >= nb) continue;  // after range

        mat(i, j) = data[lo + shift];
      }
    }
  }

  return mat;
}

/******************************************************************************/

// [[Rcpp::export]]
NumericMatrix access_dense_subset_corr_compact(Environment X,
                                               const IntegerVector& ind_row,
                                               const IntegerVector& ind_col) {

  XPtr<SFBM_corr_compact> sfbm = X["address"];
  const NumericVector p = X["p"];
  const IntegerVector first_i = X["first_i"];
  const int16_t * data = sfbm->i_x();
  const IntegerVector ind_row2 = ind_row - 1;

  int n = ind_row.size();
  int m = ind_col.size();
  NumericMatrix mat(n, m);

  for (int j = 0; j < m; j++) {

    int j2 = ind_col[j] - 1;

    int min_i = first_i[j2];
    if (min_i >= 0) {

      size_t lo = p[j2];
      int nb = p[j2 + 1] - lo;

      for (int i = 0; i < n; i++) {

        int i2 = ind_row2[i];
        if (i2 < min_i) continue;  // before range
        int shift = i2 - min_i;
        if (shift >= nb) continue;  // after range

        int16_t val = data[lo + shift];
        mat(i, j) = val / 32767.0;
      }
    }
  }

  return mat;
}

/******************************************************************************/
