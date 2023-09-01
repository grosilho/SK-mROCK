#include "OdeProblems.h"
#include <random>

DahlquistTestProblem::DahlquistTestProblem()
{
    problem_name = "DahlquistTestProblem";
    
    tend=0.1;
    
    neqn=1;
    cte_rho = true;
    know_rho = true;
    analytical_df=true;
    dense_Jacobian=true;
    
    lambda = -100.;
    xi = -2.;
}

Real DahlquistTestProblem::lambda;
Real DahlquistTestProblem::xi;

void DahlquistTestProblem::set_initial_value(Vector& y0)
{ 
    y0(0)=1.;
}

void DahlquistTestProblem::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = (lambda+xi)*x(0);
}

void DahlquistTestProblem::AN_df(Real t, Vector& x, Matrix& dfx)
{   
    dfx.resize(neqn,neqn);
    dfx(0,0) = lambda+xi;
}

void DahlquistTestProblem::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = abs(lambda+xi);
}


ScalarNonStiffNonLinearTest::ScalarNonStiffNonLinearTest()
{
    problem_name = "ScalarNonStiffNonLinearTest";
    
    tend=10.;
    
    neqn=1;
    cte_rho = false;
    know_rho = true;
    analytical_df=true;
    dense_Jacobian=true;
}

void ScalarNonStiffNonLinearTest::set_initial_value(Vector& y0)
{ 
    y0(0)=0.;
}

void ScalarNonStiffNonLinearTest::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = 0.25*x(0)+0.5*sqrt(x(0)*x(0)+1.);
}

void ScalarNonStiffNonLinearTest::AN_df(Real t, Vector& x, Matrix& dfx)
{   
    dfx.resize(neqn,neqn);
    dfx(0,0) = 0.25+0.5*x(0)/sqrt(x(0)*x(0)+1);
}

void ScalarNonStiffNonLinearTest::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = abs(0.25+0.5*y(0)/sqrt(y(0)*y(0)+1.))+1;
}

// Neuron cable equation

NeuronCable::NeuronCable()
{
    problem_name = "NeuronCable";
    
    tend=10.;
    nu = 0.01;
    beta= 1.0;
    
    neqn=128;
    cte_rho = true;
    know_rho = true;
    analytical_df=false;
    dense_Jacobian=false;
}

Real NeuronCable::nu;
Real NeuronCable::beta;

void NeuronCable::set_initial_value(Vector& y0)
{    
    const Real pi=4.*atan(1.);
    for (int j=0;j<neqn;j++)
    {
        Real x=((Real)j)/(neqn-1.);
        y0(j)=-70.+20.*cos(15.*pi*x)*(1.-x);
    }
}

void NeuronCable::f(Real t, Vector& x, Vector& fx)
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

void NeuronCable::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = nu*4.*(neqn-1)*(neqn-1)+beta;
}

// This is the brusseltaor in 2d
Brusselator::Brusselator()
{
    problem_name = "Brusselator";
    
    tend=40.0;
    
    neqn=2;
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian=true;
    
    alpha = 1.9;
}

Real Brusselator::alpha;

void Brusselator::set_initial_value(Vector& y0)
{    
    y0(0)= -0.1;
    y0(1)= 0;
}

void Brusselator::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = (alpha-1.)*x(0)+alpha*x(0)*x(0)+(x(0)+1.)*(x(0)+1.)*x(1);
    fx(1) = -alpha*x(0)-alpha*x(0)*x(0)-(x(0)+1.)*(x(0)+1.)*x(1);
}

PDEBrusselator::PDEBrusselator()
{
    problem_name = "PDEBrusselator";
    
    tend=15.0;
    
    Nu = 257; // Must be odd
    Nv = (Nu-1)/2;
    
    neqn = Nu+Nv;
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian=false;
    
    A = 1.;
    B = 3.;
    alpha = 1./50.;
    Hu = 1./(Nu+1.);
    Hv = 1./(Nv+1.);
}

Real PDEBrusselator::A;
Real PDEBrusselator::B;
Real PDEBrusselator::alpha;
unsigned int PDEBrusselator::Nu;
unsigned int PDEBrusselator::Nv;
Real PDEBrusselator::Hu;
Real PDEBrusselator::Hv;


void PDEBrusselator::set_initial_value(Vector& y0)
{    
    const Real pi=4.*atan(1.);
    
    for(int i=0;i<Nu;i++)            
        y0(i) = 1.+sin(2.*pi*(i+1.)*Hu);

    for(int i=Nu;i<Nu+Nv;i++)
        y0(i) = 3.;
}

void PDEBrusselator::f(Real t, Vector& x, Vector& fx)
{   
    Real ul = 1.;
    Real ur = 1.;
    Real vl = 3.;
    Real vr = 3.;
    int mod,j1,j2;
    
    fx(0) = A + x(0)*x(0)*(x(Nu)+vl)/2. - (B+1.)*x(0) + alpha*(ul-2.*x(0)+x(1))/Hu/Hu;
    fx(Nu-1) = A + x(Nu-1)*x(Nu-1)*(x(Nu+Nv-1)+vr)/2. - (B+1.)*x(Nu-1) + alpha*(x(Nu-2)-2.*x(Nu-1)+ur)/Hu/Hu;
    for(int i=1;i<Nu-1;i++)
    {
        mod = i % 2;
        j2 = Nu + (i-mod)/2;
        j1 = j2-1+mod;
        fx(i) = A + x(i)*x(i)*(x(j1)+x(j2))/2. - (B+1.)*x(i) + alpha*(x(i-1)-2.*x(i)+x(i+1))/Hu/Hu;
    }
    
    fx(Nu) = B*x(1) - x(1)*x(1)*x(Nu) + alpha*(vl-2.*x(Nu)+x(Nu+1))/Hv/Hv;
    fx(Nu+Nv-1) = B*x(Nu-2) - x(Nu-2)*x(Nu-2)*x(Nu+Nv-1) + alpha*(x(Nu+Nv-2)-2.*x(Nu+Nv-1)+vr)/Hv/Hv;
    for(int i=Nu+1;i<Nu+Nv-1;i++)
    {
        j1 = 2*(i-Nu)+1;
        fx(i) = B*x(j1) - x(j1)*x(j1)*x(i) + alpha*(x(i-1)-2.*x(i)+x(i+1))/Hv/Hv;
    }
}


PopulationDynamics::PopulationDynamics()
{
    problem_name = "PopulationDynamics";
    
    neqn = 2;
        
    cte_rho = false;
    know_rho = true;
    analytical_df=false;
    dense_Jacobian=true;
    
    tend=1.0;
    
    lambda1 = -500.;
    lambda2 = -4.;
    alpha = 1.;
}

void PopulationDynamics::set_initial_value(Vector& y0)
{    
    y0(0)=0.95;
    y0(1)=0.95;    
}

void PopulationDynamics::f(Real t, Vector& x, Vector& fx)
{   
    fx(0)=alpha*(x(1)-1.)-lambda1*x(0)*(1-x(0));
    fx(1)=-lambda2*x(1)*(1-x(1));
}

void PopulationDynamics::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax=550.;//max(abs(lambda2*(1.-2.*y(1))),abs(lambda1*(1.-2.*y(0))));
}
