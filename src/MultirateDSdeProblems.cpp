#include "MultirateDSdeProblems.h"

// MultirateDSdeDahlquistTestProblem
MultirateDSdeDahlquistTestProblem::MultirateDSdeDahlquistTestProblem()
:DSdeDahlquistTestProblem(), MultirateDSde()
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