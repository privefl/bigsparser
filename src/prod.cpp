/******************************************************************************/

#include <bigsparser/SFBM.h>
#include <bigsparser/SFBM-corr-compact.h>

using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
NumericVector prodVec(Environment X, const NumericVector& y) {
  XPtr<SFBM> sfbm = X["address"];
  return sfbm->prod(y);
}

// [[Rcpp::export]]
NumericVector cprodVec(Environment X, const NumericVector& y) {
  XPtr<SFBM> sfbm = X["address"];
  return sfbm->cprod(y);
}

/******************************************************************************/

// [[Rcpp::export]]
NumericVector corr_prodVec(Environment X, const NumericVector& y) {
  XPtr<SFBM_corr_compact> sfbm = X["address"];
  return sfbm->prod(y);
}

// [[Rcpp::export]]
NumericVector corr_cprodVec(Environment X, const NumericVector& y) {
  XPtr<SFBM_corr_compact> sfbm = X["address"];
  return sfbm->cprod(y);
}

/******************************************************************************/
