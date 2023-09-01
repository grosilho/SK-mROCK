#ifndef LEGENDREMETHODS_H
#define LEGENDREMETHODS_H

#include "MainHeader.h"


namespace LegendreMethods
{
    void CoefficientsRKL1(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa, 
                          vector<Real>& c, vector<Real>& d, unsigned int s, Real eps=0.0);
    void CoefficientsRKL2(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa,
                          vector<Real>& gamma, vector<Real>& c, unsigned int s, Real eps=0.0);
    
    Real ls_RKL1(unsigned int s, Real eps=0.0);
    Real ls_RKL2(unsigned int s, Real eps=0.0);
    
    Real L(Real x, unsigned int s);
    Real dL(Real x, unsigned int s);
    Real ddL(Real x, unsigned int s);
};

#endif /* LEGENDREMETHODS_H */

