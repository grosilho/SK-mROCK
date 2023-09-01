#ifndef NBEMPP_CLASS_FSERIES_H_
#define NBEMPP_CLASS_FSERIES_H_

#include <cmath>
#include <iostream>
#include <string>

#include "Eigen/Dense"
#include "fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

namespace uBidomain {
/**
 *  \brief wraps Fourier series and interpolations/evaluations into a common
 *         interface.
 *         An even number of coefficients is assumed in the representations
 *         [a_M,...,a_0,b_1,...,b_{M-1}], where 2 * M = N is the number of
 *         coefficients. The corresponding trigonometric polynomial reads
 *         a_0/2+\sum_{k=1}^{M-1} a_k\cos(2\pi k x) + b_k\sin(2\pi k x)
 *         + a_M/2 cos(2\pi Mx). However, a_M=0 is enforced in view of other
 *         aspects of the code
 **/
class FSeries {
 public:
  void set_FCoeff(const Eigen::VectorXd &FCoeff);
  void eval(int n_pts = -1);
  Eigen::VectorXd eval(const Eigen::VectorXd &angle) const;
  double eval_one_pt(const double t0) const;
  FSeries &interp(const Eigen::VectorXd &FVal);
  FSeries &truncate(double tol, int Mmin);
  FSeries &diff();
  const Eigen::VectorXd &get_FCoeff(void) const;
  const Eigen::VectorXd &get_FVal(void) const;
  

 private:
  Eigen::VectorXd FCoeff_;
  Eigen::VectorXd FVal_;
};
}  // namespace uBidomain
#endif
