#ifndef MULTIRATESDEPROBLEMS_H
#define	MULTIRATESDEPROBLEMS_H

#include "DSdeProblems.h"

class MultirateDSdeDahlquistTestProblem: public DSdeDahlquistTestProblem, public MultirateDSde
{
public:
    MultirateDSdeDahlquistTestProblem();
    virtual ~MultirateDSdeDahlquistTestProblem(){};

    void fF(Real t, Vector& x, Vector& fx);
    void fS(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmaxF, Real& eigmaxS);
};

#endif	/* MULTIRATESDEPROBLEMS_H */