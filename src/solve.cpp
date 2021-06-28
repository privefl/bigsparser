/******************************************************************************/

#include <bigsparser/EigenMatrixReplacement.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/******************************************************************************/

template <class C>
Rcpp::NumericVector sp_solve_sym_eigen0(const MatrixReplacement<C>& A,
                                        const Eigen::VectorXd& b,
                                        double tol,
                                        int maxiter) {

  // Solve Ax = b using iterative solvers with matrix-free version
  // Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper,
  //                          Eigen::IdentityPreconditioner> solver;
  Eigen::MINRES<MatrixReplacement<C>, Eigen::Lower | Eigen::Upper,
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

// [[Rcpp::export]]
Rcpp::NumericVector sp_solve_sym_eigen(Rcpp::Environment X,
                                       const Eigen::VectorXd& b,
                                       const Eigen::VectorXd& add_to_diag,
                                       double tol,
                                       int maxiter) {

  if (X.exists("first_i")) {
    Rcpp::XPtr<SFBM_compact> sfbm = X["address"];
    MatrixReplacement<SFBM_compact> A(sfbm, add_to_diag);
    return sp_solve_sym_eigen0(A, b, tol, maxiter);
  } else {
    Rcpp::XPtr<SFBM> sfbm = X["address"];
    MatrixReplacement<SFBM> A(sfbm, add_to_diag);
    return sp_solve_sym_eigen0(A, b, tol, maxiter);
  }
}

/******************************************************************************/
