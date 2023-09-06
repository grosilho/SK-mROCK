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

MultirateManyDiffusionTerms::MultirateManyDiffusionTerms()
:ManyDiffusionTerms(), MultirateDSde()
{
    FR.resize(Wsize);
    FR<<0,0,0,0,1,0,1,0,1,0;
    SR=Eigen::VectorXd::Ones(Wsize)-FR;
    cF = FR.cwiseProduct(c);
    cS = SR.cwiseProduct(c);
}

MultirateManyDiffusionTerms::~MultirateManyDiffusionTerms()
{
}

void MultirateManyDiffusionTerms::fF(Real t, Vector& X, Vector& fx)
{
    R << X(1),X(2),X(2),X(4),X(0)*X(3),X(1),X(4)*(X(4)-1)/2,X(5),X(0)*X(5),X(6);
    fx = Nu*(cF.cwiseProduct(R));
}

void MultirateManyDiffusionTerms::fS(Real t, Vector& X, Vector& fx)
{
    R << X(1),X(2),X(2),X(4),X(0)*X(3),X(1),X(4)*(X(4)-1)/2,X(5),X(0)*X(5),X(6);
    fx = Nu*(cS.cwiseProduct(R));
}
    