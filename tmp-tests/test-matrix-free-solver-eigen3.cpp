/******************************************************************************/

#include <bigsparser/EigenMatrixReplacement.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/******************************************************************************/

#define SOLVE(Solver) {                                                        \
Solver<MatrixReplacement, Eigen::Lower | Eigen::Upper,                         \
       Eigen::IdentityPreconditioner> solver;                                  \
solver.setMaxIterations(10 * sfbm->ncol());                                    \
solver.setTolerance(tol);                                                      \
solver.compute(A);                                                             \
x = solver.solve(b);                                                           \
std::cout << "#iterations: " << solver.iterations() << std::endl;              \
std::cout << "Estimated error: " << solver.error() << std::endl;               \
}                                                                              \

/******************************************************************************/

// [[Rcpp::export]]
Rcpp::NumericVector spsolve(Rcpp::Environment X, const Eigen::VectorXd& b,
                            double tol = 1e-10, bool use_CG = true) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  MatrixReplacement A(sfbm);

  Eigen::VectorXd x;

  // Solve Ax = b using iterative solvers with matrix-free version
  if (use_CG) {
    SOLVE(Eigen::ConjugateGradient)
  } else {
    SOLVE(Eigen::MINRES)
  }

  return Rcpp::wrap(x);
}

/******************************************************************************/
