/******************************************************************************/

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
List access_subset(Rcpp::Environment X,
                   const IntegerVector& ind_row,
                   const IntegerVector& ind_col) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  const IntegerVector ind_row2 = ind_row - 1;
  int N = sfbm->nrow();
  int n = ind_row.size();
  int m = ind_col.size();
  const double * data = sfbm->i_x();

  std::vector<int>    new_i;
  std::vector<double> new_x;
  NumericVector new_p(m + 1);

  NumericVector cur_col(N);

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

    // second pass: read requested indices, and store values (if != 0)
    for (int i = 0; i < n; i++) {
      double val = cur_col[ind_row2[i]];
      if (val != 0) {
        new_i.push_back(i);
        new_x.push_back(val);
      }
    }

    // third pass: reset column
    for (size_t k = lo; k < up; k += 2) {
      int ind = data[k];
      cur_col[ind] = 0;
    }

    // store cumulative number of elements
    new_p[j + 1] = new_i.size();
  }

  return List::create(_["p"] = new_p, _["i"] = new_i, _["x"] = new_x);
}

/******************************************************************************/

// [[Rcpp::export]]
List access_subset_compact(Rcpp::Environment X,
                           const IntegerVector& ind_row,
                           const IntegerVector& ind_col) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  const NumericVector p = X["p"];
  const IntegerVector first_i = X["first_i"];
  const IntegerVector ind_row2 = ind_row - 1;
  int n = ind_row.size();
  int m = ind_col.size();
  const double * data = sfbm->i_x();

  std::vector<int>    new_i;
  std::vector<double> new_x;
  NumericVector new_p(m + 1);

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

        double val = data[lo + shift];
        if (val != 0) {
          new_i.push_back(i);
          new_x.push_back(val);
        }
      }
    }

    // store cumulative number of elements
    new_p[j + 1] = new_i.size();
  }

  return List::create(_["p"] = new_p, _["i"] = new_i, _["x"] = new_x);
}

/******************************************************************************/
