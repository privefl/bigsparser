/******************************************************************************/

#include <bigsparser/SFBM.h>

using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/******************************************************************************/

// [[Rcpp::export]]
NumericVector prodVec(Environment X, const NumericVector& y) {

  if (X.exists("first_i")) {
    XPtr<SFBM_compact> sfbm = X["address"];
    return sfbm->prod(y);
  } else {
    XPtr<SFBM> sfbm = X["address"];
    return sfbm->prod(y);
  }
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector cprodVec(Environment X, const NumericVector& y) {

  if (X.exists("first_i")) {
    XPtr<SFBM_compact> sfbm = X["address"];
    return sfbm->cprod(y);
  } else {
    XPtr<SFBM> sfbm = X["address"];
    return sfbm->cprod(y);
  }
}

/******************************************************************************/
