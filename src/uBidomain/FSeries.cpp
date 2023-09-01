#include "FSeries.h"

namespace uBidomain {

    void FSeries::set_FCoeff(const Eigen::VectorXd &FCoeff) {
        eigen_assert(FCoeff(0) == 0 && "a_M must be zero in this implementation");
        if (FCoeff.size() % 2) {
            std::cerr << "FSeries: number of Fourier coefficients is odd\n";
            return;
        }
        FCoeff_ = FCoeff;
        return;
    }

    void FSeries::eval(int n_pts) {
        if (!(FCoeff_.size())) {
            std::cerr << "FSeries: initialize Fourier coefficients first\n";
            return;
        }
        if (n_pts == -1) n_pts = FCoeff_.size();
        FVal_.resize(n_pts);
        FVal_.setZero();
        ifft((double *) FVal_.data(), (double *) FCoeff_.data(), FCoeff_.size(), n_pts);
        return;
    }

    // angle should be in [0,1]!

    Eigen::VectorXd FSeries::eval(const Eigen::VectorXd &angle) const {
        int M = 0;
        Eigen::VectorXd FCsin;
        Eigen::VectorXd FCcos;
        Eigen::VectorXd retVal;
        double s1 = 0;
        double s2 = 0;        

        if (!(FCoeff_.size())) {
            std::cerr << "FSeries: initialize Fourier coefficients first\n";
            // FIXME return some error to stop execution
        }
        // set up Fourier coefficients in adequate format for Reinsch's algorithm
        M = FCoeff_.size() / 2;
        FCcos.resize(M + 1);
        FCsin.resize(M);
        FCcos = FCoeff_.head(M + 1).reverse();
        FCsin.tail(M - 1) = FCoeff_.tail(M - 1);
        FCcos(0) *= 0.5;
        FCcos(M) *= 0.5;
        FCsin(0) = 0;
        // evaluate Fourier seriese
        retVal.resize(angle.size());
        for (auto i = 0; i < angle.size(); ++i) {
            ReinschEval(&s1, &s2, (double *) FCsin.data(), FCsin.size(), angle(i));
            retVal(i) = s2;
            ReinschEval(&s1, &s2, (double *) FCcos.data(), FCcos.size(), angle(i));
            retVal(i) += s1;
        }
        return retVal;
    }
    
    double FSeries::eval_one_pt(const double t0) const
    {       
        Eigen::VectorXd angle(1);        
        angle(0)=t0;
//        while(angle(0)<0.)
//            angle(0)+=1.;
//        while(angle(0)>1.)
//            angle(0)-=1.;
        return eval(angle)(0);
    }

    FSeries& FSeries::interp(const Eigen::VectorXd &FVal) {
        if (FVal.size() % 2) {
            std::cerr << "FSeries: number of interpolation points is odd\n";
            return *this;
        }
        FVal_ = FVal;
        FCoeff_.resize(FVal_.size());
        FCoeff_.setZero();
        fft((double *) FVal_.data(), (double *) FCoeff_.data(), FVal_.size());
        FCoeff_(0) = 0;
        return *this;
    }

    FSeries& FSeries::truncate(double tol, int Mmin) {
        double norm2 = FCoeff_.squaredNorm();
        double cur_norm2 = 0;
        int trunck = 0;
        int M = FCoeff_.size() / 2;
        cur_norm2 = FCoeff_(M) * FCoeff_(M);
        while (trunck < M - 1 && sqrt(norm2 - cur_norm2) / sqrt(norm2) > tol) {
            ++trunck;
            cur_norm2 += FCoeff_(M - trunck) * FCoeff_(M - trunck) +
                    FCoeff_(M + trunck) * FCoeff_(M + trunck);
        }
        if (trunck == M - 1 || Mmin>M)
            return *this;
        else {
            ++trunck;
            trunck = std::max(trunck, Mmin);
            FCoeff_ = FCoeff_.segment(M - trunck, 2 * trunck);
            FCoeff_(0) = 0;
        }
        return *this;
    }

    FSeries& FSeries::diff() {
        int M = FCoeff_.size() / 2;
        // note that differentiation kills a_M in order to avoid increasing the
        // number of terms every time this thing is differentiated
        FCoeff_(0) = 0;
        FCoeff_(M) = 0;
        for (auto i = 1; i < FCoeff_.size(); ++i) FCoeff_(i) *= 2. * M_PI * (i - M);

        FCoeff_.tail(FCoeff_.size() - 1).reverseInPlace();
        return *this;
    }

    const Eigen::VectorXd& FSeries::get_FCoeff(void) const {
        return FCoeff_;
    }

    const Eigen::VectorXd& FSeries::get_FVal(void) const {
        return FVal_;
    }

} // namespace NBEM
