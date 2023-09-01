#include <cmath>
#include "DSdeProblems.h"

DSdeDahlquistTestProblem::DSdeDahlquistTestProblem()
:Ode(), DSde(), DahlquistTestProblem()       
{
    Wsize = 1;
    noise=DIAGONAL;
    
    sigma  = 10.;//0.7*sqrt(2.*abs(lambda+xi));
}

Real DSdeDahlquistTestProblem::phi(const Vector& X)
{
    return sin(X(0));
}

void DSdeDahlquistTestProblem::g(Real t, Vector& x, Vector& G)
{
    G(0) = sigma*x(0);
}

void DSdeDahlquistTestProblem::g(Real t, Vector& x, Vector& G, int r)
{
    G(0) = sigma*x(0);
}

void DSdeDahlquistTestProblem::g(Real t, Vector& x, Matrix& G)
{
    G(0,0) = sigma*x(0);
}

// -----------------------------------------------------------------------------

DSdeScalarNonStiffNonLinearTest::DSdeScalarNonStiffNonLinearTest()
:Ode(), DSde(), ScalarNonStiffNonLinearTest()       
{
    Wsize = 1;
    noise=DIAGONAL;
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
    
    neqn=1;
    Wsize = 10;
    noise=GENERAL;
    
    cte_rho = true;
    know_rho = true;
    
    a.resize(Wsize);
    b.resize(Wsize);
    a<< 10., 15., 20., 25., 40., 25., 20., 15., 20., 25.;
    b<< 2., 4., 5., 10., 20., 2., 4., 5., 10., 20.;
    
    for(int i=0;i<Wsize;i++)
    {
        a(i)=1./a(i);
        b(i)=1./b(i);
    }
    
}

ManyDiffusionTerms::~ManyDiffusionTerms()
{
}

void ManyDiffusionTerms::set_initial_value(Vector& y0)
{    
    y0(0)=1.;
}

void ManyDiffusionTerms::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = x(0);
}

void ManyDiffusionTerms::analytical_df(Real t, Vector& x, Matrix& dfx)
{
    dfx.resize(neqn,neqn);
    dfx(0,0) = 1.;
}

void ManyDiffusionTerms::g(Real t, Vector& x, Vector& G, int r)
{   
    G(0) = a(r-1)*sqrt(x(0)+b(r-1));
}

void ManyDiffusionTerms::g(Real t, Vector& x, Matrix& G)
{
    for(int r=0;r<Wsize;r++)
        G(0,r) = a(r)*sqrt(x(0)+b(r));
}

Real ManyDiffusionTerms::phi(const Vector& X)
{
    return X(0)*X(0);
}

void ManyDiffusionTerms::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = 1.;
}

/*
SdeOneDimensionalTest3::SdeOneDimensionalTest3(string output_file)
:Ode(string("Tests/SdeOneDimensionalTest3/")+output_file),
 Sde(string("Tests/SdeOneDimensionalTest3/")+output_file)
{
    neqn = 1;
    Wsize = 1;
    
    noise=DIAGONAL;
    
    cte_rho = false;
    know_rho = true;
    cte_df=false;
    stiff=false;
    
    y0=0.;
    t0=0.0;
    tend=1.0;
    t=t0;
    
    y = new Vector(neqn);
    
    init_solution();
//    output_solution();
}

void SdeOneDimensionalTest3::init_solution()
{    
    (*y)(0)=y0;
}

void SdeOneDimensionalTest3::f(Real t, Vector& x, Vector& fx)
{   
    Real co = cos(x(0));
    fx(0) = -sin(x(0))*co*co*co;
}

void SdeOneDimensionalTest3::g(Real t, Vector& x, Vector& gx)
{
    Real co = cos(x(0));
    gx(0) = co*co;
}

Real SdeOneDimensionalTest3::phi()
{
    return tan((*y)(0))*tan((*y)(0));
}

Real SdeOneDimensionalTest3::Exact_phi()
{
    return tan(y0)*tan(y0)+tend;
}

void SdeOneDimensionalTest3::rho(Real& eigmax)
{
    Real co = cos((*y)(0));
    eigmax = abs( co*co*(2.*cos(2*(*y)(0))-1.));
}

SdeOneDimensionalTest4::SdeOneDimensionalTest4(string output_file)
:Ode(string("Tests/SdeOneDimensionalTest4/")+output_file),
 Sde(string("Tests/SdeOneDimensionalTest4/")+output_file)
{
    neqn = 1;
    Wsize = 1;
    
    noise=DIAGONAL;
    
    cte_rho = false;
    know_rho = true;
    cte_df=false;
    stiff=true;
    
    t0=0.0;
    tend=1.0;
    t=t0;
    alpha = 0.5;
    
    y = new Vector(neqn);
    
    init_solution();
//    output_solution();
}

void SdeOneDimensionalTest4::init_solution()
{    
    (*y)(0)=1.;
}

void SdeOneDimensionalTest4::f(Real t, Vector& x, Vector& fx)
{   
    Real eps = 1e-10;
    if(t<alpha)
        fx(0)=0.;
    else
        fx(0) = 0.5*x(0)/sqrt(t-alpha+eps);
}

void SdeOneDimensionalTest4::g(Real t, Vector& x, Vector& gx)
{
    gx(0) = x(0);
}

Real SdeOneDimensionalTest4::phi()
{
    return (*y)(0);
}

Real SdeOneDimensionalTest4::Exact_phi()
{
    return exp(sqrt(tend-alpha));
}

void SdeOneDimensionalTest4::rho(Real& eigmax)
{
    Real eps = 1e-10;
    if(t<alpha)
        eigmax=1.;
    else
        eigmax = 0.5/sqrt(t-alpha+eps);
}

TwoDimensionalSinCos::TwoDimensionalSinCos(string output_file)
:Ode(string("Tests/TwoDimensionalSinCos/")+output_file),
 Sde(string("Tests/TwoDimensionalSinCos/")+output_file)
{
    neqn = 2;
    Wsize = 2;
    
    noise=GENERAL;
    
    cte_rho = true;
    know_rho = true;
    cte_df=true;
    stiff=false;
    
    t0=0.0;
    tend=1.0;
    t=t0;
    
    y = new Vector(neqn);
    
    init_solution();
//    output_solution();
}

void TwoDimensionalSinCos::init_solution()
{    
    (*y)(0)=1.;
    (*y)(1)=1.;
}

void TwoDimensionalSinCos::f(Real t, Vector& x, Vector& fx)
{   
    fx(0)=-x(1);
    fx(1)=x(0);
}

void TwoDimensionalSinCos::g(Real t, Vector& x, Matrix& G)
{
    G(0,0) = 0.;
    G(1,0) = sin(x(0)+x(1))/sqrt(1.+t);
    
    G(0,1) = cos(x(0)+x(1))/sqrt(1.+t);
    G(1,1) = 0.;
}

void TwoDimensionalSinCos::g(Real t, Vector& x, Vector& G, int r)
{
    if(r==0)
    {
        G(0) = 0.;
        G(1) = sin(x(0)+x(1))/sqrt(1.+t);
    }
    else if(r==1)
    {
        G(0) = cos(x(0)+x(1))/sqrt(1.+t);
        G(1) = 0.;
    }
    else
        cout<<"Error in computing diffusion"<<endl;
}

Real TwoDimensionalSinCos::phi()
{
    return (*y)(0)*(*y)(0)+(*y)(1)*(*y)(1);
}

Real TwoDimensionalSinCos::Exact_phi()
{
    return 2.+log(1.+tend);
}

void TwoDimensionalSinCos::rho(Real& eigmax)
{
    eigmax= 1.;
}


DuffingVanDerPol::DuffingVanDerPol(string output_file)
:Ode(string("Tests/DuffingVanDerPol/")+output_file),
 Sde(string("Tests/DuffingVanDerPol/")+output_file)
{
    neqn = 2;
    Wsize = 1;
    
    noise=GENERAL;
    
    cte_rho = false;
    know_rho = true;
    cte_df=false;
    stiff=false;
    
    t0=0.0;
    tend=8.0;
    t=t0;
    
    alpha = 1.;
    sigma = 1.;
    
    y = new Vector(neqn);
    
    init_solution();
//    output_solution();
}

void DuffingVanDerPol::init_solution()
{    
    (*y)(0)=-3.;
    (*y)(1)=0.;
}

void DuffingVanDerPol::f(Real t, Vector& x, Vector& fx)
{   
    fx(0)= x(1);
    fx(1)= x(0)*(alpha-x(0)*x(0))-x(1);
}

void DuffingVanDerPol::g(Real t, Vector& x, Matrix& G)
{
    G(0,0) = 0.;
    G(1,0) = sigma*x(0);
}

void DuffingVanDerPol::g(Real t, Vector& x, Vector& G, int r)
{
    if(r==0)
    {
        G(0) = 0.;
        G(1) = sigma*x(0);
    }
    else
        cout<<"Error in computing diffusion"<<endl;
}

Real DuffingVanDerPol::phi()
{
    return (*y)(0)*(*y)(0)+(*y)(1)*(*y)(1);
}

Real DuffingVanDerPol::Exact_phi()
{
    //reference solution
    //1e9 iterations, dt=1e-2
    return 1.00997637240578;
}

void DuffingVanDerPol::rho(Real& eigmax)
{
    if(4.*alpha-12*(*y)(0)*(*y)(0)+1.>=0)
        eigmax = 0.5*(1.+sqrt(4.*alpha-12*(*y)(0)*(*y)(0)+1.));
    else
        eigmax = 0.5*sqrt(-4*alpha+12*(*y)(0)*(*y)(0));
}

StochasticBrusselator::StochasticBrusselator(string output_file)
:Ode(string("Tests/StochasticBrusselator/")+output_file),
 Sde(string("Tests/StochasticBrusselator/")+output_file),
 Brusselator(string("Tests/StochasticBrusselator/")+output_file)       
{
    Wsize = 1;
    noise=GENERAL;
    
    sigma = 0.1;
    
    //init solution in other class
}

Real StochasticBrusselator::phi()
{
    return (*y)(0)+(*y)(1);
}

Real StochasticBrusselator::Exact_phi()
{
    return 0.;
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


StochasticPopulationDynamics::StochasticPopulationDynamics(string output_file)
:Ode(string("Tests/StochasticPopulationDynamics/")+output_file),
 Sde(string("Tests/StochasticPopulationDynamics/")+output_file),
 PopulationDynamics(string("Tests/StochasticPopulationDynamics/")+output_file)
{
    Wsize = 1;
    
    noise=GENERAL;
    
    mu1 = sqrt(abs(lambda1));
    mu2 = 1.;
}

void StochasticPopulationDynamics::g(Real t, Vector& x, Matrix& G)
{
    G(0,0) = -mu1*x(0)*(1.-x(0));
    G(1,0) = -mu2*x(1)*(1.-x(1));
}

void StochasticPopulationDynamics::g(Real t, Vector& x, Vector& gx, int r)
{
    gx(0) = -mu1*x(0)*(1.-x(0));
    gx(1) = -mu2*x(1)*(1.-x(1));
}

Real StochasticPopulationDynamics::phi()
{
    return (*y)(0)+(*y)(1);
}

Real StochasticPopulationDynamics::Exact_phi()
{
    return 1.9990206435918350;
}


StochasticNeuronCable::StochasticNeuronCable(string output_file)
:Ode(string("Tests/StochasticNeuronCable/")+output_file),
 Sde(string("Tests/StochasticNeuronCable/")+output_file),
 NeuronCable(string("Tests/StochasticNeuronCable/")+output_file)       
{
    Wsize = 128;
    noise=DIAGONAL;
    
    sigma = 4e-3;
    V0 = 10.;
    
    //init solution in other class
}

Real StochasticNeuronCable::phi()
{
    return y->sum();
//    return accu(*y);
}

Real StochasticNeuronCable::Exact_phi()
{
    return -3.050088673263552e+03;
}

void StochasticNeuronCable::g(Real t, Vector& x, Vector& G)
{
    Real sqrtdx=1./sqrt(neqn-1.);
    
    for(int i=0;i<neqn;i++)
        G(i)=sigma*(x(i)+V0)/sqrtdx;
}
 */