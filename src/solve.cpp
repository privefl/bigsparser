/******************************************************************************/

#include <bigsparser/EigenMatrixReplacement.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/******************************************************************************/

// [[Rcpp::export]]
Rcpp::NumericVector sp_solve_sym_eigen(Rcpp::Environment X,
                                       const Eigen::VectorXd& b,
                                       const Eigen::VectorXd& add_to_diag,
                                       double tol,
                                       int maxiter) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  MatrixReplacement A(sfbm, add_to_diag);

  // Solve Ax = b using iterative solvers with matrix-free version
  // Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper,
  //                          Eigen::IdentityPreconditioner> solver;
  Eigen::MINRES<MatrixReplacement, Eigen::Lower | Eigen::Upper,
                Eigen::IdentityPreconditioner> solver;

  solver.setTolerance(tol);
  solver.setMaxIterations(maxiter);

  solver.compute(A);
  Eigen::VectorXd x = solver.solve(b);

  double eps = solver.error();
  // std::cout << eps << " in " << solver.iterations() << " iterations." << std::endl;

  if (std::isnan(eps))
    Rcpp::stop("Solver failed.");
  if (eps > tol)
    Rcpp::warning("Estimated error: %s.", eps);

  return Rcpp::wrap(x);
}

/******************************************************************************/
