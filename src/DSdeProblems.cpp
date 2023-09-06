#include <cmath>
#include "DSdeProblems.h"

DSdeDahlquistTestProblem::DSdeDahlquistTestProblem()
:Ode(), DSde()
{
    problem_name = "DSdeDahlquistTestProblem";
    
    tend=0.1;
    
    neqn=1;
    cte_rho = true;
    know_rho = true;
    analytical_df=true;
    dense_Jacobian=true;
    
    lambda = -100.;
    mu = -10;
    xi  = 0.7*sqrt(2.*abs(lambda+mu));

    Wsize = 1;
    noise=DIAGONAL;    
}

Real DSdeDahlquistTestProblem::lambda;
Real DSdeDahlquistTestProblem::mu;
Real DSdeDahlquistTestProblem::xi;

void DSdeDahlquistTestProblem::set_initial_value(Vector& y0)
{ 
    y0(0)=1.;
}

void DSdeDahlquistTestProblem::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = (lambda+mu)*x(0);
}

void DSdeDahlquistTestProblem::AN_df(Real t, Vector& x, Matrix& dfx)
{   
    dfx.resize(neqn,neqn);
    dfx(0,0) = lambda+mu;
}

void DSdeDahlquistTestProblem::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = abs(lambda+mu);
}

Real DSdeDahlquistTestProblem::phi(const Vector& X)
{
    return sin(X(0));
}

void DSdeDahlquistTestProblem::g(Real t, Vector& x, Vector& G)
{
    G(0) = xi*x(0);
}

void DSdeDahlquistTestProblem::g(Real t, Vector& x, Vector& G, int r)
{
    G(0) = xi*x(0);
}

void DSdeDahlquistTestProblem::g(Real t, Vector& x, Matrix& G)
{
    G(0,0) = xi*x(0);
}

// -----------------------------------------------------------------------------

DSdeScalarNonStiffNonLinearTest::DSdeScalarNonStiffNonLinearTest()
:Ode(), DSde()
{
    problem_name = "DSdeScalarNonStiffNonLinearTest";
    
    tend=10.;
    
    neqn=1;
    cte_rho = false;
    know_rho = true;
    analytical_df=true;
    dense_Jacobian=true;

    Wsize = 1;
    noise=DIAGONAL;
}

void DSdeScalarNonStiffNonLinearTest::set_initial_value(Vector& y0)
{ 
    y0(0)=0.;
}

void DSdeScalarNonStiffNonLinearTest::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = 0.25*x(0) + 0.5*sqrt(x(0)*x(0)+1.);
}

void DSdeScalarNonStiffNonLinearTest::AN_df(Real t, Vector& x, Matrix& dfx)
{   
    dfx.resize(neqn,neqn);
    dfx(0,0) = 0.25+0.5*x(0)/sqrt(x(0)*x(0)+1);
}

void DSdeScalarNonStiffNonLinearTest::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = abs(0.25+0.5*y(0)/sqrt(y(0)*y(0)+1));
}

Real DSdeScalarNonStiffNonLinearTest::phi(const Vector& X)
{
    return asinh(X(0));
}

void DSdeScalarNonStiffNonLinearTest::g(Real t, Vector& x, Vector& G)
{
    G(0) = sqrt((x(0)*x(0)+1.)/2.);
}

void DSdeScalarNonStiffNonLinearTest::g(Real t, Vector& x, Vector& G, int r)
{
    G(0) = sqrt((x(0)*x(0)+1.)/2.);
}

void DSdeScalarNonStiffNonLinearTest::g(Real t, Vector& x, Matrix& G)
{
    G(0,0) = sqrt((x(0)*x(0)+1.)/2.);
}


ManyDiffusionTerms::ManyDiffusionTerms()
:Ode(), DSde()
{
    problem_name = "ManyDiffusionTerms";
    
    tend=1.;
    
    neqn=7;
    Wsize = 10;
    noise = GENERAL;
    
    cte_rho = false;
    know_rho = false;
    
    R.resize(Wsize);
    c.resize(Wsize);
    c<< 0.0078, 0.043, 0.0039, 0.0007, 0.38, 3, 0.5, 5, 0.12, 9;

    Nu.resize(neqn,Wsize);
    Nu << 1.,0.,0.,0.,-1.,1.,0.,0.,-1.,1.,
        -1.,0.,0.,0.,1.,-1.,0.,0.,0.,0.,
         1.,0.,-1.,0.,0.,0.,0.,0.,0.,0.,
         1.,0.,0.,0.,-1.,1.,0.,0.,0.,0.,
         0.,1.,0.,-1.,0.,0.,-2.,2.,0.,0.,
         0.,0.,0.,0.,0.,0,1.,-1.,-1.,1.,
         0.,0.,0.,0.,0.,0.,0.,0.,1.,-1.;    
}

ManyDiffusionTerms::~ManyDiffusionTerms()
{
}

void ManyDiffusionTerms::set_initial_value(Vector& y0)
{    
    y0.resize(neqn);
    y0<<0.,1.,0.,30.,0.,0.,0.;
    y0 *= 5000.;
}

void ManyDiffusionTerms::f(Real t, Vector& X, Vector& fx)
{   
    R << X(1),X(2),X(2),X(4),X(0)*X(3),X(1),X(4)*(X(4)-1)/2,X(5),X(0)*X(5),X(6);
    fx = Nu*(c.cwiseProduct(R));
}

void ManyDiffusionTerms::g(Real t, Vector& X, Vector& G, int r)
{   
    R << X(1),X(2),X(2),X(4),X(0)*X(3),X(1),X(4)*(X(4)-1)/2,X(5),X(0)*X(5),X(6);
    G = Nu.col(r)*sqrt(abs(c(r)*R(r)));
}

void ManyDiffusionTerms::g(Real t, Vector& X, Matrix& G)
{
    R << X(1),X(2),X(2),X(4),X(0)*X(3),X(1),X(4)*(X(4)-1)/2,X(5),X(0)*X(5),X(6);
    D.diagonal() << c.cwiseProduct(R).cwiseAbs().cwiseSqrt();
    G = Nu*D;
}

Real ManyDiffusionTerms::phi(const Vector& X)
{
    return X(0)*X(0);
}


StochasticBrusselator::StochasticBrusselator()
: Ode(), DSde()
{
    problem_name = "StochasticBrusselator";
    
    tend=40.0;
    
    neqn=2;
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian=true;

    Wsize = 1;
    noise=GENERAL;
    
    alpha = 1.9;
    sigma = 1.;
}

void StochasticBrusselator::set_initial_value(Vector& y0)
{    
    y0(0)= -0.1;
    y0(1)= 0;
}

void StochasticBrusselator::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = (alpha-1.)*x(0)+alpha*x(0)*x(0)+(x(0)+1.)*(x(0)+1.)*x(1);
    fx(1) = -alpha*x(0)-alpha*x(0)*x(0)-(x(0)+1.)*(x(0)+1.)*x(1);
}

Real StochasticBrusselator::phi(const Vector& X)
{
    return X(0)+X(1);
}

void StochasticBrusselator::g(Real t, Vector& x, Vector& G, int r)
{
    G(0) = sigma*x(0)*(1.+x(0));
    G(1) = -sigma*x(0)*(1.+x(0));
    
    if(r>0)
        cout<<"Problem in computing diffusion"<<endl;
}

void StochasticBrusselator::g(Real t, Vector& x, Matrix& G)
{
    G(0,0) = sigma*x(0)*(1.+x(0));
    G(1,0) = -sigma*x(0)*(1.+x(0));
}

StochasticNeuronCable::StochasticNeuronCable()
:Ode(), DSde()      
{
    problem_name = "StochasticNeuronCable";
    
    tend=10.;
    nu = 0.01;
    beta= 1.0;
    sigma = 0.5;
    V0 = 10.;
    
    neqn=128;
    cte_rho = true;
    know_rho = true;
    analytical_df=false;
    dense_Jacobian=false;

    Wsize = 128;
    noise=DIAGONAL;
}

StochasticNeuronCable::~StochasticNeuronCable()
{
}

void StochasticNeuronCable::set_initial_value(Vector& y0)
{    
    const Real pi=4.*atan(1.);
    for (int j=0;j<neqn;j++)
    {
        Real x=((Real)j)/(neqn-1.);
        y0(j)=-70.+20.*cos(15.*pi*x)*(1.-x);
    }
}

void StochasticNeuronCable::f(Real t, Vector& x, Vector& fx)
{   
    // Computing diffusion with Neumann bnd conditions
    fx(0)=nu*2.*(x(1)-x(0))*(neqn-1.)*(neqn-1.)-beta*x(0);
    fx(neqn-1)=nu*2.*(x(neqn-2)-x(neqn-1))*(neqn-1.)*(neqn-1.)-beta*x(neqn-1);
    for (int i=1;i<neqn-1;i++)
    {
        fx(i)=nu*(x(i-1)-2.*x(i)+x(i+1))*(neqn-1.)*(neqn-1.)-beta*x(i);
        if(abs(i/(neqn-1.)-0.5)<0.1)
            fx(i) += 5.*exp(1.-1e-2/(1e-2-(i/(neqn-1.)-0-5)*(i/(neqn-1.)-0.5)));
    }
}

void StochasticNeuronCable::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = nu*4.*(neqn-1)*(neqn-1)+beta;
}

Real StochasticNeuronCable::phi(const Vector& X)
{
    return X.sum();
}

void StochasticNeuronCable::g(Real t, Vector& x, Vector& G)
{
    Real sqrtdx=1./sqrt(neqn-1.);
    
    for(int i=0;i<neqn;i++)
        G(i)=sigma*(x(i)+V0)/sqrtdx;
}
