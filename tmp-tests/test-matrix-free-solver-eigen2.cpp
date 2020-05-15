// https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(bigsparser)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#define STRICT_R_HEADERS

#include <RcppEigen.h>
#include <bigsparser/SFBM.h>

// #include <iostream>
// #include <Eigen/Core>
// #include <Eigen/Dense>
// #include <Eigen/IterativeLinearSolvers>
// #include <unsupported/Eigen/IterativeSolvers>

class MatrixReplacement;

namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template<>
struct traits<MatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
{};
}
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  Index rows() const { return sfbm->nrow(); }
  Index cols() const { return sfbm->ncol(); }

  template<typename Rhs>
  Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // Custom API:
  MatrixReplacement(SFBM * sfbm) : sfbm(sfbm) {}

  SFBM * matrix() const { return sfbm; }

private:
  SFBM * sfbm;
};


// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {

template<typename Rhs>
struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixReplacement,Rhs,generic_product_impl<MatrixReplacement,Rhs> >
{
  typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;

  template<typename Dest>
  static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
  {
    // This method should implement "dst += alpha * lhs * rhs" inplace,
    // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
    assert(alpha==Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    dst.noalias() += (lhs.matrix())->prod<Eigen::VectorXd>(rhs);
  }
};

}
}

// [[Rcpp::export]]
Rcpp::NumericVector spsolve(Rcpp::Environment X, const Eigen::VectorXd& b, double tol = 1e-10) {

  Rcpp::XPtr<SFBM> sfbm = X["address"];
  MatrixReplacement A(sfbm);

  Eigen::VectorXd x;

  // Solve Ax = b using various iterative solver with matrix-free version:
  {
    Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
    cg.setMaxIterations(10 * sfbm->ncol());
    cg.setTolerance(tol);
    cg.compute(A);
    x = cg.solve(b);
    std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    // std::cout << x << std::endl;
  }

  {
    Eigen::BiCGSTAB<MatrixReplacement, Eigen::IdentityPreconditioner> bicg;
    bicg.setMaxIterations(10 * sfbm->ncol());
    bicg.setTolerance(tol);
    bicg.compute(A);
    x = bicg.solve(b);
    std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
    // std::cout << x << std::endl;
  }

  {
    Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
    gmres.setMaxIterations(10 * sfbm->ncol());
    gmres.setTolerance(tol);
    gmres.compute(A);
    x = gmres.solve(b);
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    // std::cout << x << std::endl;
  }

  {
    Eigen::DGMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
    gmres.setMaxIterations(10 * sfbm->ncol());
    gmres.setTolerance(tol);
    gmres.compute(A);
    x = gmres.solve(b);
    std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    // std::cout << x << std::endl;
  }

  {
    Eigen::MINRES<MatrixReplacement, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
    minres.setMaxIterations(10 * sfbm->ncol());
    minres.setTolerance(tol);
    minres.compute(A);
    x = minres.solve(b);
    std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;
    // std::cout << x << std::endl;
  }

  return Rcpp::wrap(x);
}
