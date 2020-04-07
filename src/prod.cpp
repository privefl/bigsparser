/******************************************************************************/

#include <bigsparser/SFBM.h>

using namespace Rcpp;

/******************************************************************************/

template<class C>
C SFBM::prod(const C& x) {

  C res(n);
  for (int i = 0; i < n; i++) res[i] = 0;

  for (int j = 0; j < m; j++) {

    auto lo = data + p[j];
    auto up = data + p[j + 1];

    for (auto it = lo; it != up; it++) {
      res[it->i] += it->x * x[j];
    }
  }

  return res;
}

template<class C>
C SFBM::cprod(const C& x) {

  C res(m);

  for (int j = 0; j < m; j++)
    res[j] = dot_col(j, x);

  return res;
}

template<class C>
double SFBM::dot_col(int j, const C& x) {

  auto lo = data + p[j];
  auto up = data + p[j + 1];

  double cp = 0;
  for (auto it = lo; it != up; it++)
    cp += it->x * x[it->i];

  return cp;
}

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
