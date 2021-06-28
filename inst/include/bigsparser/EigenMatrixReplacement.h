#ifndef EIGEN_MATRIX_REPLACEMENT_H
#define EIGEN_MATRIX_REPLACEMENT_H

/******************************************************************************/

// https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(bigsparser)]]

#define STRICT_R_HEADERS

#include <RcppEigen.h>
#include <bigsparser/SFBM.h>

template <class C> class MatrixReplacement;

/******************************************************************************/

namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template <class C> struct traits< MatrixReplacement<C> > :
  public Eigen::internal::traits< Eigen::SparseMatrix<double> > {};
}
}

/******************************************************************************/

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
template <class C>
class MatrixReplacement : public Eigen::EigenBase< MatrixReplacement<C> > {
public:
  // Required typedefs, constants, and method:
  typedef typename Eigen::EigenBase<MatrixReplacement<C> >::Index Index;
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
  Eigen::Product<MatrixReplacement<C>, Rhs, Eigen::AliasFreeProduct>
  operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixReplacement<C>, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // Custom API:
  MatrixReplacement(C * sfbm) :
    sfbm(sfbm), add_to_diag(Eigen::VectorXd::Zero(sfbm->ncol())) {}

  MatrixReplacement(C * sfbm, const Eigen::VectorXd& add_to_diag) :
    sfbm(sfbm), add_to_diag(add_to_diag) {}

  C * matrix() const { return sfbm; }
  Eigen::VectorXd extra() const { return add_to_diag; }

private:
  C * sfbm;
  const Eigen::VectorXd add_to_diag;
};

/******************************************************************************/

// Implementation of MatrixReplacement * Eigen::DenseVector
// though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {

template<class C, typename Rhs>
struct generic_product_impl<MatrixReplacement<C>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixReplacement<C>,Rhs,generic_product_impl<MatrixReplacement<C>,Rhs> >
{
  typedef typename Product<MatrixReplacement<C>, Rhs>::Scalar Scalar;

  template<typename Dest>
  static void scaleAndAddTo(Dest& dst, const MatrixReplacement<C> & lhs, const Rhs& rhs, const Scalar& alpha)
  {
    // This method should implement "dst += alpha * lhs * rhs" inplace,
    // however, for iterative solvers, alpha is always equal to 1,
    // so let's not bother about it.
    assert(alpha==Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    // Use cprod because matrix is symmetric (and should be faster than prod)
    // template -> https://stackoverflow.com/a/37995805/6103040
    C * sfbm_ptr = lhs.matrix();
    dst.noalias() += sfbm_ptr->template cprod<Eigen::VectorXd>(rhs) +
      lhs.extra().cwiseProduct(rhs);
  }
};

}
}

/******************************************************************************/

#endif // EIGEN_MATRIX_REPLACEMENT_H
