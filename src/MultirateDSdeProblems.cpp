#include "MultirateDSdeProblems.h"

MultirateDSdeDahlquistTestProblem::MultirateDSdeDahlquistTestProblem()
:DSdeDahlquistTestProblem(), MultirateDSde()
{
}

MultirateDSdeDahlquistTestProblem::~MultirateDSdeDahlquistTestProblem()
{
}

void MultirateDSdeDahlquistTestProblem::fF(Real t, Vector& x, Vector& fx)
{
    fx(0) = lambda*x(0);
}

void MultirateDSdeDahlquistTestProblem::fS(Real t, Vector& x, Vector& fx)
{
    fx(0) = mu*x(0);
}
    
void MultirateDSdeDahlquistTestProblem::rho(Real t, Vector& y, Real& eigmaxF, Real& eigmaxS)
{
    eigmaxF = abs(lambda);
    eigmaxS = abs(mu);
}

MultirateDSdeScalarNonStiffNonLinearTest::MultirateDSdeScalarNonStiffNonLinearTest()
:DSdeScalarNonStiffNonLinearTest(), MultirateDSde()
{
}

MultirateDSdeScalarNonStiffNonLinearTest::~MultirateDSdeScalarNonStiffNonLinearTest()
{
}

void MultirateDSdeScalarNonStiffNonLinearTest::fF(Real t, Vector& x, Vector& fx)
{
    fx(0) = 0.5*sqrt(x(0)*x(0)+1.);
}

void MultirateDSdeScalarNonStiffNonLinearTest::fS(Real t, Vector& x, Vector& fx)
{
    fx(0) = 0.25*x(0);
}
    
void MultirateDSdeScalarNonStiffNonLinearTest::rho(Real t, Vector& y, Real& eigmaxF, Real& eigmaxS)
{
    eigmaxF = abs(0.5*y(0)/sqrt(y(0)*y(0)+1));
    eigmaxS = 0.25;
}