#ifndef EIGEN_MATRIX_REPLACEMENT_H
#define EIGEN_MATRIX_REPLACEMENT_H

/******************************************************************************/

// https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(bigsparser)]]

#define STRICT_R_HEADERS

#include <RcppEigen.h>
#include <bigsparser/SFBM.h>

class MatrixReplacement;

/******************************************************************************/

namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template<> struct traits<MatrixReplacement> :
  public Eigen::internal::traits< Eigen::SparseMatrix<double> > {};
}
}

/******************************************************************************/

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
  Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>
  operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // Custom API:
  MatrixReplacement(SFBM * sfbm) :
    sfbm(sfbm), add_to_diag(Eigen::VectorXd::Zero(sfbm->ncol())) {}

  MatrixReplacement(SFBM * sfbm, const Eigen::VectorXd& add_to_diag) :
    sfbm(sfbm), add_to_diag(add_to_diag) {}

  SFBM * matrix() const { return sfbm; }
  Eigen::VectorXd extra() const { return add_to_diag; }

private:
  SFBM * sfbm;
  const Eigen::VectorXd add_to_diag;
};

/******************************************************************************/

// Implementation of MatrixReplacement * Eigen::DenseVector
// though a specialization of internal::generic_product_impl:
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
    // however, for iterative solvers, alpha is always equal to 1,
    // so let's not bother about it.
    assert(alpha==Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    dst.noalias() += (lhs.matrix())->prod<Eigen::VectorXd>(rhs) +
      (lhs.extra()).cwiseProduct(rhs);
  }
};

}
}

/******************************************************************************/

#endif // EIGEN_MATRIX_REPLACEMENT_H
