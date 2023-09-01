#ifndef MATRIXREPLACEMENT_H
#define MATRIXREPLACEMENT_H

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "Ode.h"
 
class MatrixReplacement;
using Eigen::SparseMatrix;
 
namespace Eigen {
namespace internal {
  // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<MatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};
}
}
 
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
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
 
  Index rows() const { return M; }
  Index cols() const { return N; }
 
  template<typename Rhs>
  Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
 
  // Custom API:
  MatrixReplacement(Ode* ode_, Real t_, Real h_, 
                    Vector* x_, Vector* tmp1_, Vector* tmp2_);
  const Real get_t() const;
  const Real get_h() const;
  Ode* get_ode() const;
  Vector* get_x() const;
  Vector* get_tmp1() const;
  Vector* get_tmp2() const;
  
  Ode* ode;
  unsigned int M,N;
  const Real h;
  const Real t;
  Vector* x;
  Vector* tmp1;
  Vector* tmp2;
          
};

MatrixReplacement::MatrixReplacement(Ode* ode_, Real t_, Real h_, 
        Vector* x_, Vector* tmp1_, Vector* tmp2_)
:ode(ode_), t(t_), h(h_), x(x_), tmp1(tmp1_), tmp2(tmp2_)
{
    M = ode->get_system_size();
    N = ode->get_system_size();
}

const Real MatrixReplacement::get_h() const
{
    return h;
}

const Real MatrixReplacement::get_t() const
{
    return t;
}

Ode* MatrixReplacement::get_ode() const
{
    return ode;
}

Vector* MatrixReplacement::get_x() const
{
    return x;
}

Vector* MatrixReplacement::get_tmp1() const
{
    return tmp1;
}

Vector* MatrixReplacement::get_tmp2() const
{
    return tmp2;
}
 
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
      // Here I implement the matrix multilication by I-h*J, with J the Jacobian of ode->f
      
        Vector* x = lhs.get_x();
        Vector* tmp1 = lhs.get_tmp1();
        Vector* tmp2 = lhs.get_tmp2();
        Real eps = 1e-8*x->norm()/rhs.norm();
        Real h = lhs.get_h();
        Real t = lhs.get_t();
      
        lhs.get_ode()->f(t+h,*x,*tmp1);
        lhs.get_ode()->f(t+h,*x+eps*rhs,*tmp2);
        
        dst -= alpha*h*(*tmp2-*tmp1)/eps;
        dst += alpha*rhs;
      
    }
    
  };
 
}
}

#endif /* MATRIXREPLACEMENT_H */

