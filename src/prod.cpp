/******************************************************************************/

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector prodVec(Environment X, const NumericVector& y) {
  XPtr<SFBM> sfbm = X["address"];
  return sfbm->prod(y);
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector cprodVec(Environment X, const NumericVector& y) {
  XPtr<SFBM> sfbm = X["address"];
  return sfbm->cprod(y);
}

/******************************************************************************/
