#ifndef CHEBYSHEVMETHODS_H
#define CHEBYSHEVMETHODS_H

#include "MainHeader.h"


namespace ChebyshevMethods
{
    void CoefficientsRKC1(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa, 
                          vector<Real>& c, vector<Real>& d, unsigned int s, Real eps=0.05);
    void CoefficientsRKC2(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa,
                          vector<Real>& gamma, vector<Real>& c, unsigned int s, Real eps=0.05);
    void CoefficientsSKROCK(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa, 
                          vector<Real>& c, unsigned int s, Real eps=0.05);
    
    void CoefficientsRKU1(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa, 
                          vector<Real>& c, vector<Real>& d, unsigned int s, Real eps=0.05);
    void CoefficientsRKU2(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa,
                          vector<Real>& gamma, vector<Real>& c, unsigned int s, Real eps=0.05);
    
    Real ls_RKC1(unsigned int s, Real eps=0.0);
    Real ls_RKC2(unsigned int s, Real eps=0.0);
    Real ls_SKROCK(unsigned int s, Real eps=0.0);
    Real ls_RKU1(unsigned int s, Real eps=0.0);
    Real ls_RKU2(unsigned int s, Real eps=0.0);
    
    Real T(Real x, unsigned int s);
    Real dT(Real x, unsigned int s);
    Real ddT(Real x, unsigned int s);
    Real U(Real x, unsigned int s);
    Real dU(Real x, unsigned int s);
    Real ddU(Real x, unsigned int s);
};

#endif /* CHEBYSHEVMETHODS_H */

