#ifndef MULTIRATESDEPROBLEMS_H
#define	MULTIRATESDEPROBLEMS_H

#include "DSdeProblems.h"

class MultirateDSdeDahlquistTestProblem: public DSdeDahlquistTestProblem, public MultirateDSde
{
public:
    MultirateDSdeDahlquistTestProblem();
    virtual ~MultirateDSdeDahlquistTestProblem();

    void fF(Real t, Vector& x, Vector& fx);
    void fS(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmaxF, Real& eigmaxS);
};

class MultirateDSdeScalarNonStiffNonLinearTest: public DSdeScalarNonStiffNonLinearTest, public MultirateDSde
{
public:
    MultirateDSdeScalarNonStiffNonLinearTest();
    virtual ~MultirateDSdeScalarNonStiffNonLinearTest();

    void fF(Real t, Vector& x, Vector& fx);
    void fS(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmaxF, Real& eigmaxS);
};

class MultirateManyDiffusionTerms: public ManyDiffusionTerms, public MultirateDSde
{
public:
    MultirateManyDiffusionTerms();
    virtual ~MultirateManyDiffusionTerms();
    
    void fF(Real t, Vector& x, Vector& fx);
    void fS(Real t, Vector& x, Vector& fx);
        
protected:
    Vector cF,cS;
    Vector FR,SR;
};

#endif	/* MULTIRATESDEPROBLEMS_H */