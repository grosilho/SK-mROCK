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
    autodiff_df = true;
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
    autodiff_df = true;
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
    autodiff_df=true;
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
    autodiff_df=false;
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
    autodiff_df = true;
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


// Test problem 10 from Krogh
Krogh10::Krogh10()
{
    problem_name = "Krogh10";
    
    tend = 6.19216933131963970674;
    mu = 1./82.45;
    mus = 1.-mu;
    
    neqn=4;
    cte_rho = false;
    know_rho = true;
    analytical_df=false;
    autodiff_df = true;
    dense_Jacobian=true;
}

Real Krogh10::mu;
Real Krogh10::mus;

void Krogh10::set_initial_value(Vector& y0)
{    
    // These initial conditions give a periodic solution, so they are also
    // exact solution at tend
    y0(0)=1.2;
    y0(1)=0.0;
    y0(2)=0.0;
    y0(3)= -1.0493575098031990726;
}

void Krogh10::f(Real t, Vector& x, Vector& fx)
{   
    Real y1 = x(0);
    Real y2 = x(1);
    Real y1p = x(2);
    Real y2p = x(3);
    Real r1 = sqrt(((y1+mu)*(y1+mu)+y2*y2));
    Real r2 = sqrt(((y1-mus)*(y1-mus)+y2*y2));
    r1 = r1*r1*r1;
    r2 = r2*r2*r2;
    
    fx(0) = y1p;
    fx(1) = y2p;
    fx(2) = 2.*y2p+y1-mus*(y1+mu)/r1-mu*(y1-mus)/r2;
    fx(3) = -2.*y1p+y2-mus*y2/r1-mu*y2/r2;
}

void Krogh10::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = abs(0.25+0.5*y(0)/sqrt(y(0)*y(0)+1.));
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

VanDerPol::VanDerPol()
{
    problem_name = "VanDerPol";
    
    tend = 2.0;
    
    eps=1e-4;
    
    neqn=2;
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian=true;
}

void VanDerPol::set_initial_value(Vector& y0)
{    
    y0(0)=2.;
    y0(1)=-0.6;
}

void VanDerPol::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = x(1);
    fx(1) = ((1.-x(0)*x(0))*x(1)-x(0))/eps;
}

ODEIonicModel::ODEIonicModel()
{
    problem_name = "ODEIonicModel";
    
//    im = std::make_unique<uBidomain::HodgkinHuxley>();
//    im = std::make_unique<uBidomain::MitchellSchaeffer>();
//    im = std::make_unique<uBidomain::Cubic>();
    im = std::make_unique<uBidomain::MyoKit::MyoKitIonicModel>();
    
    tend = im->get_tend();
    
    n_gating_vars = im->get_n_gating_vars();
    neqn = 1+n_gating_vars;
    
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian=true;
}

void ODEIonicModel::set_initial_value(Vector& y0)
{    
    im->get_y0(y0);
}

void ODEIonicModel::f(Real t, Vector& y, Vector& fy)
{           
    fy = im->f(t,y);
}

DiffusionRefinedMesh::DiffusionRefinedMesh()
{
    problem_name = "DiffusionRefinedMesh";
    
    tend=1.;
    nu = 0.1;
    
    N1 = 100; // N1 points in [0,1]
    N2 = 200; // N2 points in ]1,2]    
    
    if(N2<N1)
    {
        cout<<"Errore N2<N1"<<endl;
        return;
    }
    
    neqn= N1+N2+1;
    cte_rho = true;
    know_rho = true;
    analytical_df=false;
    dense_Jacobian = false;
    
    H1 = 1./N1;
    H2 = 1./N2;
}

void DiffusionRefinedMesh::set_initial_value(Vector& y0)
{    
    const Real pi=4.*atan(1.);
    Real x;
    for (int j=0;j<N1;j++)
    {
        x = j*H1;
        y0(j)=sin(x*pi);
    }
    for (int j=N1;j<=N1+N2;j++)
    {
        x = 1 + (j-N1)*H2;
        y0(j)=sin(x*pi);
    }
}

void DiffusionRefinedMesh::f(Real t, Vector& x, Vector& fx)
{   
    fx(0)=nu*2.*(x(1)-x(0))/H1/H1; // Neumann bnd cond
    fx(neqn-1)=nu*2.*(x(neqn-2)-x(neqn-1))/H2/H2;
    fx(N1) = nu*(H2*x(N1-1)-(H1+H2)*x(N1)+H1*x(N1+1))/(0.5*H1*H2*(H1+H2));
    for(int i=1;i<N1;i++)
        fx(i)=nu*(x(i-1)-2.*x(i)+x(i+1))/H1/H1;
    for(int i=N1+1;i<N1+N2;i++)
        fx(i)=nu*(x(i-1)-2.*x(i)+x(i+1))/H2/H2;
}

void DiffusionRefinedMesh::rho(Real t, Vector& y, Real& eigmax)
{
    eigmax = nu*4.*max(1./H1/H1,1./H2/H2);
}

InfectiousDiseaseTransmission::InfectiousDiseaseTransmission()
{
    problem_name = "InfectiousDiseaseTransmission";
    
    tend = 50.;
       
    neqn=3;
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian=true;
    
    o = 1e4;
    e = 1e4;
    m = 0.01;
    p = 2e5;
    g = 50.;
    r = 0.8;
    ID50 = 1e5; //1e5
    a = log2(exp(1))/ID50;
}

void InfectiousDiseaseTransmission::set_initial_value(Vector& y0)
{    
    //y = (S,I,V)
    y0(0) = 5e6; //1e7
    y0(1) = 1e6; //1e6
    y0(2) = ID50;
}

void InfectiousDiseaseTransmission::f(Real t, Vector& z, Vector& fz)
{           
    Real S,I,V, fV;
    S = z(0);
    I = z(1);
    V = z(2);
    
    fV = 1.-exp(-a*V);
    
    fz(0) = p-m*S-r*S*fV;
    fz(1) = r*S*fV-(m+g)*I;
    fz(2) = o*I-e*V;
    
    if(abs(t-25.)<5.)
        fz(2) += 1e5;
}

RobertsonChemicalSystem::RobertsonChemicalSystem()
{
    problem_name = "RobertsonChemicalSystem";
    
    tend = 1e2;
       
    neqn=3;
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    autodiff_df=true;
    dense_Jacobian=true;
    
    k1 = 0.04;
    k2 = 1e4;
    k3 = 3e7;
}

Real RobertsonChemicalSystem::k1;
Real RobertsonChemicalSystem::k2;
Real RobertsonChemicalSystem::k3;

void RobertsonChemicalSystem::set_initial_value(Vector& y0)
{    
    y0(0) = 1.;
    y0(1) = 2e-5;
    y0(2) = 0.1;
}

void RobertsonChemicalSystem::f(Real t, Vector& z, Vector& fz)
{                   
    fz(0) = -k1*z(0) + k2*z(1)*z(2);
    fz(1) =  k1*z(0) - k2*z(1)*z(2)-k3*z(1)*z(1);
    fz(2) =                         k3*z(1)*z(1);
}

Oregonator::Oregonator()
{
    problem_name = "Oregonator";
    
    tend=360.0;
    
    neqn=5;
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian=true;
    
    k1 = 1.34;
    k2 = 1.6e9;
    k3 = 8e3;
    k4 = 4e7;
    k5 = 1.;
    ff=1.;
    A = 0.06;
    B = A;
}

void Oregonator::set_initial_value(Vector& y0)
{   
    // y = (Br03-, Br-, HBrO2, P, Ce(IV))
    y0(0)= 0.6e-1;
    y0(1)= 0.33e-6;
    y0(2)= 0.501e-10;
    y0(3)= 0.3e-1;
    y0(4)= 0.24e-7;
}

void Oregonator::f(Real t, Vector& x, Vector& fx)
{   
    fx(0) = -k1*x(0)*x(1) - k3*x(0)*x(2);
    fx(1) = -k1*x(0)*x(1) - k2*x(1)*x(2) + k5*x(4);
    fx(2) =  k1*x(0)*x(1) - k2*x(1)*x(2) + k3*x(0)*x(2) - 2.*k4*x(2)*x(2);   
    fx(3) =  k2*x(1)*x(2) + k4*x(2)*x(2);
    fx(4) =  k3*x(0)*x(2) - k5*x(4);
}

CUSP::CUSP()
{
    problem_name = "CUSP";
    
    tend=2.0;
    
    neqn=3*128;
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian = false;
    
    eps = 1e-4;
    sigma = 1./144.;
}

void CUSP::set_initial_value(Vector& y0)
{    
    int n=neqn/3;
    const Real pi=4.*atan(1.);
    
    for(int i=0;i<n;i++)
    {
        Real x = ((Real)i)/n;
        y0(i) = 0.;
        y0(i+n) = -2.*cos(2.*pi*x);
        y0(i+2*n) = 2*sin(2.*pi*x);
    }
}

void CUSP::f(Real t, Vector& x, Vector& fx)
{   
    int n=neqn/3;
    Real D = sigma*n*n;
    
    fx(0) = D*(x(n-1)-2.*x(0)+x(1)) - (pow(x(0),3)+x(n)*x(0)+x(2*n))/eps ;
    fx(n-1) = D*(x(n-2)-2.*x(n-1)+x(0)) - (pow(x(n-1),3)+x(2*n-1)*x(n-1)+x(3*n-1))/eps;
    for(int i=1;i<n-1;i++)
        fx(i) = D*(x(i-1)-2.*x(i)+x(i+1)) - (pow(x(i),3)+x(n+i)*x(i)+x(2*n+i))/eps;
    
    fx(n) = D*(x(2*n-1)-2.*x(n)+x(n+1)) + (x(2*n)+0.07*v(x(0)));
    fx(2*n-1) = D*(x(2*n-2)-2.*x(2*n-1)+x(n)) + (x(3*n-1)+0.07*v(x(n-1)));
    for(int i=n+1;i<2*n-1;i++)
        fx(i) = D*(x(i-1)-2.*x(i)+x(i+1)) + (x(i+n)+0.07*v(x(i-n)));
    
    fx(2*n) = D*(x(3*n-1)-2.*x(2*n)+x(2*n+1)) + (1.-pow(x(n),2))*x(2*n)-x(n)-0.4*x(0)+0.035*v(x(0));
    fx(3*n-1) = D*(x(3*n-2)-2.*x(3*n-1)+x(2*n)) + (1.-pow(x(2*n-1),2))*x(3*n-1)-x(2*n-1)-0.4*x(n-1)+0.035*v(x(n-1));
    for(int i=2*n+1;i<3*n-1;i++)
        fx(i) = D*(x(i-1)-2.*x(i)+x(i+1)) + (1.-pow(x(i-n),2))*x(i)-x(i-n)-0.4*x(i-2*n)+0.035*v(x(i-2*n));
}

Real CUSP::v(Real y)
{
    Real u = (y-0.7)*(y-1.3);
    return u/(u+0.1);
}


mReactionDiffusion2DEquations::mReactionDiffusion2DEquations()
{
    problem_name = "mReactionDiffusion2DEquations";
    
    tend=1.0;
    
    N = 100; // Mesh NxN
    m = 4;   // m reactants -> m Reaction-Diffusion equations on NxN mesh
    Nsq = N*N;
    
    neqn = m*Nsq;
    
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian = false;
    
    D.resize(m);
    if(m==4)
    {
        //Fig 2a
    //    D[0] = 16.7;
    //    D[1] = 49.5;
    //    D[2] = 36.4;
    //    D[3] = 117.6;
        //Fig 3a
    //    D[0] = 12.6;
    //    D[1] = 47.5;
    //    D[2] = 27.5;
    //    D[3] = 141.5;
    //    D[4] = 1.;
        //Fig 3b
//        D[0] = 1.85;
//        D[1] = 50.6;
//        D[2] = 5.66;
//        D[3] = 186.;
        //Fig 3c
        D[0] = 1.31;
        D[1] = 34.;
        D[2] = 9.87;
        D[3] = 3449; //originally 344.9

        alpha = 1.; // 1 for Fig 3bc, 0.1 else
        beta = alpha;
        A = 3.; //always 3
        B = 6.; //6 for Fig 3c, 9 else
    }
    
    h = 1./N;
}

void mReactionDiffusion2DEquations::set_initial_value(Vector& y0)
{    

//    for(unsigned int k=0;k<m;k++)
//        for(unsigned int i=0;i<N;i++)
//            for(unsigned int j=0;j<N;j++)
//                y(k*Nsq+i*N+j) = sin((k+1)*2*3.141592*(i*i+j*j)*h*h); // sin((k+1)2pi (x^1+y^2))     
    
    default_random_engine generator;
    uniform_real_distribution<double> distribution(-0.2,0.2);
  
    if(m==4)
    {
        for(unsigned int k=0;k<2;k++)
        for(unsigned int i=0;i<N;i++)
        for(unsigned int j=0;j<N;j++)
            y0(k*Nsq + N*i + j) = A + distribution(generator);
        
        for(unsigned int k=2;k<4;k++)
        for(unsigned int i=0;i<N;i++)
        for(unsigned int j=0;j<N;j++)
            y0(k*Nsq + N*i + j) = B/A + distribution(generator);
    }

}

void mReactionDiffusion2DEquations::f(Real t, Vector& y, Vector& fy)
{   
    unsigned int kNsq;
    
    // Intra-layer diffusion
    for(unsigned int k=0; k<m; k++)
    {
        kNsq = k*Nsq;
        fy(kNsq) =  D[k]*(y(kNsq+N-1)-2.*y(kNsq)+y(kNsq+1))/h/h                                     //Lower-Left corner
                  + D[k]*(y(kNsq+N)-2.*y(kNsq)+y(kNsq+Nsq-N))/h/h;
        fy(kNsq+N-1) =  D[k]*(y(kNsq)-2.*y(kNsq+N-1)+y(kNsq+N-2))/h/h                               //Lower-Right corner
                      + D[k]*(y(kNsq+2*N-1)-2.*y(kNsq+N-1)+y(kNsq+Nsq-1))/h/h;
        fy(kNsq+Nsq-N) =  D[k]*(y(kNsq+Nsq-1)-2.*y(kNsq+Nsq-N)+y(kNsq+Nsq-N+1))/h/h                 //Upper-Left corner
                        + D[k]*(y(kNsq+Nsq-2*N)-2.*y(kNsq+Nsq-N)+y(kNsq))/h/h;
        fy(kNsq+Nsq-1) =  D[k]*(y(kNsq+Nsq-2)-2.*y(kNsq+Nsq-1)+y(kNsq+Nsq-N))/h/h                   //Upper-Right corner
                        + D[k]*(y(kNsq+Nsq-N-1)-2.*y(kNsq+Nsq-1)+y(kNsq+N-1))/h/h;
        for(unsigned int j=1;j<N-1;j++)
        {
            fy(kNsq+j) =  D[k]*(y(kNsq+j-1)-2.*y(kNsq+j)+y(kNsq+j+1))/h/h                           // Lower boundary
                        + D[k]*(y(kNsq+j+N)-2.*y(kNsq+j)+y(kNsq+j+Nsq-N))/h/h;
            fy(kNsq+Nsq-N+j) =  D[k]*(y(kNsq+Nsq-N+j-1)-2.*y(kNsq+Nsq-N+j)+y(kNsq+Nsq-N+j+1))/h/h   // Upper boundary
                              + D[k]*(y(kNsq+Nsq-N+j-N)-2.*y(kNsq+Nsq-N+j)+y(kNsq+j))/h/h;
            fy(kNsq+j*N) =  D[k]*(y(kNsq+j*N+1)-2.*y(kNsq+j*N)+y(kNsq+j*N+N-1))/h/h                 // Left boundary
                          + D[k]*(y(kNsq+j*N-N)-2.*y(kNsq+j*N)+y(kNsq+j*N+N))/h/h;
            fy(kNsq+j*N+N-1) =  D[k]*(y(kNsq+j*N+N-2)-2.*y(kNsq+j*N+N-1)+y(kNsq+j*N))/h/h           // Right boundary
                              + D[k]*(y(kNsq+j*N+2*N-1)-2.*y(kNsq+j*N+N-1)+y(kNsq+j*N-1))/h/h;
            
            for(unsigned int i=1;i<N-1;i++)
                fy(kNsq+j*N+i) =  D[k]*(y(kNsq+j*N+i+1)-2.*y(kNsq+j*N+i)+y(kNsq+j*N+i-1))/h/h       // Internal points
                                + D[k]*(y(kNsq+j*N-N+i)-2.*y(kNsq+j*N+i)+y(kNsq+j*N+N+i))/h/h;
        }        
   }
    
    // Inter-layer diffusion and reaction terms
    for(unsigned int k=0; k<m; k++)    
    for(unsigned int i=0;i<N;i++)
    for(unsigned int j=0;j<N;j++)
        fy(k*Nsq + N*i + j) += inter_diff_and_reac(y, k, i,j);
    
}

Real mReactionDiffusion2DEquations::inter_diff_and_reac(Vector& y, unsigned int k, unsigned int i, unsigned int j)
{

    if(m == 4) 
    {
        unsigned int k1,k2;
        Real xk, xk2, yk, yk2;
    
        if(k<=1)
        {            
            k2 = 1-k;
            xk = y(k*Nsq+i*N+j);
            xk2 = y(k2*Nsq+i*N+j);
            yk = y((k+2)*Nsq+i*N+j);
            return alpha*(xk2-xk) + A - (1.+B)*xk + xk*xk*yk;
        }
        if(k>=2)
        {
            k2 = 5-k;
            yk = y(k*Nsq+i*N+j);
            yk2 = y(k2*Nsq+i*N+j);
            xk = y((k-2)*Nsq+i*N+j);
            return beta*(yk2-yk) + B*xk - xk*xk*yk;
        }
    }
    else
    {
        cout<<"Error in mReactionDiffusion2DEquations::inter_diff_and_reac"<<endl;
        return 0.;
    }
    return 0.;
}

RadiationDiffusion::RadiationDiffusion()
{
    problem_name = "RadiationDiffusion";
    
    tend=1.2; // originally 3.0;
    
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian = false;
    
    N = 128; // Mesh NxN
    Nsq = N*N;
    neqn = 2*Nsq;
    
    h = 1./(N-1);
    
    E0 = 1e-5;
    T0 = sqrt(sqrt(E0));
    
    k = 0.005;
}

void RadiationDiffusion::set_initial_value(Vector& y0)
{
    //index i*N+j corresponds to E(x) with x=(j*h,i*h)
    //index Nsq+i*N+j corresponds to T(x) with x=(j*h,i*h)
    for(unsigned int i=0;i<N;i++)
    for(unsigned int j=0;j<N;j++)
    {
        y0(i*N+j)=E0;
        y0(Nsq+i*N+j)=T0;
    }
}

void RadiationDiffusion::f(Real t, Vector& x, Vector& fx)
{
    //index i*N+j corresponds to E(x) with x=(j*h,i*h)
    //index Nsq+i*N+j corresponds to T(x) with x=(j*h,i*h)
    Real Eij,Eijp1,Eijm1,Eip1j,Eim1j,Eip1jp1,Eim1jp1,Eip1jm1,Eim1jm1;
    Real Tij,Tijp1,Tijm1,Tip1j,Tim1j;
    Real Eijp05,Eijm05,Eip05j,Eim05j;
    Real Tijp05,Tijm05,Tip05j,Tim05j;
    Real gradEijp05,gradEijm05,gradEip05j,gradEim05j;
    Real radiation;
    
    //Diffusion and radiation in the interior of the domain
    for(unsigned int i=1;i<N-1;i++)
    for(unsigned int j=1;j<N-1;j++)
    {
        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
        Eij = x(N*i+j);
        Eijp1 = x(N*i+j+1);
        Eijm1 = x(N*i+j-1);
        Eip1j = x(N*(i+1)+j);
        Eim1j = x(N*(i-1)+j);
        Eip1jp1 = x(N*(i+1)+j+1);
        Eim1jp1 = x(N*(i-1)+j+1);
        Eip1jm1 = x(N*(i+1)+j-1);
        Eim1jm1 = x(N*(i-1)+j-1);
        Tij = x(Nsq+N*i+j);
        Tijp1 = x(Nsq+N*i+j+1);
        Tijm1 = x(Nsq+N*i+j-1);
        Tip1j = x(Nsq+N*(i+1)+j);
        Tim1j = x(Nsq+N*(i-1)+j);
        Eijp05 = (Eij+Eijp1)/2.;
        Eijm05 = (Eij+Eijm1)/2.;
        Eip05j = (Eij+Eip1j)/2.;
        Eim05j = (Eij+Eim1j)/2.;
        Tijp05 = (Tij+Tijp1)/2.;
        Tijm05 = (Tij+Tijm1)/2.;
        Tip05j = (Tij+Tip1j)/2.;
        Tim05j = (Tij+Tim1j)/2.;
        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
        gradEijp05 = sqrt(gradEijp05);
        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
        gradEijm05 = sqrt(gradEijm05);
        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
        gradEip05j = sqrt(gradEip05j);
        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
        gradEim05j = sqrt(gradEim05j);
        
        // computing diffusion
        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                                -D2(Tijm05)*(Tij  -Tijm1)
                                +D2(Tip05j)*(Tip1j-Tij)
                                -D2(Tim05j)*(Tij  -Tim1j));
        
        //adding radiation
        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
        fx(i*N+j) += radiation;
        fx(Nsq+i*N+j) -= radiation;
    }
    
    // upper boundary
    for(unsigned int j=1;j<N-1;j++)
    {
        //same code as in the interior of the domain but we use homogeneous Neumann 
        // and all i+1 become i-1
        unsigned int i=N-1;
        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
        Eij = x(N*i+j);
        Eijp1 = x(N*i+j+1);
        Eijm1 = x(N*i+j-1);
        Eip1j = x(N*(i-1)+j);
        Eim1j = x(N*(i-1)+j);
        Eip1jp1 = x(N*(i-1)+j+1);
        Eim1jp1 = x(N*(i-1)+j+1);
        Eip1jm1 = x(N*(i-1)+j-1);
        Eim1jm1 = x(N*(i-1)+j-1);
        Tij = x(Nsq+N*i+j);
        Tijp1 = x(Nsq+N*i+j+1);
        Tijm1 = x(Nsq+N*i+j-1);
        Tip1j = x(Nsq+N*(i-1)+j);
        Tim1j = x(Nsq+N*(i-1)+j);
        Eijp05 = (Eij+Eijp1)/2.;
        Eijm05 = (Eij+Eijm1)/2.;
        Eip05j = (Eij+Eip1j)/2.;
        Eim05j = (Eij+Eim1j)/2.;
        Tijp05 = (Tij+Tijp1)/2.;
        Tijm05 = (Tij+Tijm1)/2.;
        Tip05j = (Tij+Tip1j)/2.;
        Tim05j = (Tij+Tim1j)/2.;
        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
        gradEijp05 = sqrt(gradEijp05);
        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
        gradEijm05 = sqrt(gradEijm05);
        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
        gradEip05j = sqrt(gradEip05j);
        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
        gradEim05j = sqrt(gradEim05j);
        
        // computing diffusion
        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                                -D2(Tijm05)*(Tij  -Tijm1)
                                +D2(Tip05j)*(Tip1j-Tij)
                                -D2(Tim05j)*(Tij  -Tim1j));
        
        //adding radiation
        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
        fx(i*N+j) += radiation;
        fx(Nsq+i*N+j) -= radiation;
    }
    
    // lower boundary
    for(unsigned int j=1;j<N-1;j++)
    {
        unsigned int i=0;
        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
        Eij = x(N*i+j);
        Eijp1 = x(N*i+j+1);
        Eijm1 = x(N*i+j-1);
        Eip1j = x(N*(i+1)+j);
        Eim1j = x(N*(i+1)+j);
        Eip1jp1 = x(N*(i+1)+j+1);
        Eim1jp1 = x(N*(i+1)+j+1);
        Eip1jm1 = x(N*(i+1)+j-1);
        Eim1jm1 = x(N*(i+1)+j-1);
        Tij = x(Nsq+N*i+j);
        Tijp1 = x(Nsq+N*i+j+1);
        Tijm1 = x(Nsq+N*i+j-1);
        Tip1j = x(Nsq+N*(i+1)+j);
        Tim1j = x(Nsq+N*(i+1)+j);
        Eijp05 = (Eij+Eijp1)/2.;
        Eijm05 = (Eij+Eijm1)/2.;
        Eip05j = (Eij+Eip1j)/2.;
        Eim05j = (Eij+Eim1j)/2.;
        Tijp05 = (Tij+Tijp1)/2.;
        Tijm05 = (Tij+Tijm1)/2.;
        Tip05j = (Tij+Tip1j)/2.;
        Tim05j = (Tij+Tim1j)/2.;
        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
        gradEijp05 = sqrt(gradEijp05);
        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
        gradEijm05 = sqrt(gradEijm05);
        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
        gradEip05j = sqrt(gradEip05j);
        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
        gradEim05j = sqrt(gradEim05j);
        
        // computing diffusion
        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                                -D2(Tijm05)*(Tij  -Tijm1)
                                +D2(Tip05j)*(Tip1j-Tij)
                                -D2(Tim05j)*(Tij  -Tim1j));
        
        //adding radiation
        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
        fx(i*N+j) += radiation;
        fx(Nsq+i*N+j) -= radiation;
    }
    
    //left boundary
    for(unsigned int i=1;i<N-1;i++)
    {
        unsigned int j=0;
        //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=1 for E
        //yielding that value for Eijm1
        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
        Eij = x(N*i+j);
        Eijp1 = x(N*i+j+1);        
        Eip1j = x(N*(i+1)+j);
        Eim1j = x(N*(i-1)+j);
        Eip1jp1 = x(N*(i+1)+j+1);
        Eim1jp1 = x(N*(i-1)+j+1);        
        Tij = x(Nsq+N*i+j);
        Tijp1 = x(Nsq+N*i+j+1);
        Tijm1 = x(Nsq+N*i+j+1);
        Tip1j = x(Nsq+N*(i+1)+j);
        Tim1j = x(Nsq+N*(i-1)+j);
        Eijm1 = Eijp1+3.*h*sigma(j*h,i*h,Tij)*(4.-Eij);
        Eip1jm1 = Eip1jp1+3.*h*sigma(j*h,(i+1)*h,Tip1j)*(4.-Eip1j);
        Eim1jm1 = Eim1jp1+3.*h*sigma(j*h,(i-1)*h,Tim1j)*(4.-Eim1j);
        Eijp05 = (Eij+Eijp1)/2.;
        Eijm05 = (Eij+Eijm1)/2.;
        Eip05j = (Eij+Eip1j)/2.;
        Eim05j = (Eij+Eim1j)/2.;
        Tijp05 = (Tij+Tijp1)/2.;
        Tijm05 = (Tij+Tijm1)/2.;
        Tip05j = (Tij+Tip1j)/2.;
        Tim05j = (Tij+Tim1j)/2.;
        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
        gradEijp05 = sqrt(gradEijp05);
        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
        gradEijm05 = sqrt(gradEijm05);
        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
        gradEip05j = sqrt(gradEip05j);
        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
        gradEim05j = sqrt(gradEim05j);
        
        // computing diffusion
        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                                -D2(Tijm05)*(Tij  -Tijm1)
                                +D2(Tip05j)*(Tip1j-Tij)
                                -D2(Tim05j)*(Tij  -Tim1j));
        
        //adding radiation
        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
        fx(i*N+j) += radiation;
        fx(Nsq+i*N+j) -= radiation;
    }
    
    //right boundary
    for(unsigned int i=1;i<N-1;i++)
    {
        //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=0 for E
        //yielding that value for Eijp1
        unsigned int j=N-1;
        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
        Eij = x(N*i+j);        
        Eijm1 = x(N*i+j-1);
        Eip1j = x(N*(i+1)+j);
        Eim1j = x(N*(i-1)+j);
        Eip1jm1 = x(N*(i+1)+j-1);
        Eim1jm1 = x(N*(i-1)+j-1);
        Tij = x(Nsq+N*i+j);
        Tijp1 = x(Nsq+N*i+j-1);
        Tijm1 = x(Nsq+N*i+j-1);
        Tip1j = x(Nsq+N*(i+1)+j);
        Tim1j = x(Nsq+N*(i-1)+j);
        Eijp1 = Eijm1-3.*h*sigma(j*h,i*h,Tij)*Eij;
        Eip1jp1 = Eip1jm1-3.*h*sigma(j*h,(i+1)*h,Tip1j)*Eip1j;
        Eim1jp1 = Eim1jm1-3.*h*sigma(j*h,(i-1)*h,Tim1j)*Eim1j;
        Eijp05 = (Eij+Eijp1)/2.;
        Eijm05 = (Eij+Eijm1)/2.;
        Eip05j = (Eij+Eip1j)/2.;
        Eim05j = (Eij+Eim1j)/2.;
        Tijp05 = (Tij+Tijp1)/2.;
        Tijm05 = (Tij+Tijm1)/2.;
        Tip05j = (Tij+Tip1j)/2.;
        Tim05j = (Tij+Tim1j)/2.;
        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
        gradEijp05 = sqrt(gradEijp05);
        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
        gradEijm05 = sqrt(gradEijm05);
        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
        gradEip05j = sqrt(gradEip05j);
        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
        gradEim05j = sqrt(gradEim05j);
        
        // computing diffusion
        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                                -D2(Tijm05)*(Tij  -Tijm1)
                                +D2(Tip05j)*(Tip1j-Tij)
                                -D2(Tim05j)*(Tij  -Tim1j));
        
        //adding radiation
        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
        fx(i*N+j) += radiation;
        fx(Nsq+i*N+j) -= radiation;
    }
    
    // upper left corner
    unsigned int i=N-1;
    unsigned int j=0;
    //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=1 for E
    //yielding that value for Eijm1
    //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
    Eij = x(N*i+j);
    Eijp1 = x(N*i+j+1);        
    Eip1j = x(N*(i-1)+j);
    Eim1j = x(N*(i-1)+j);
    Eip1jp1 = x(N*(i-1)+j+1);
    Eim1jp1 = x(N*(i-1)+j+1);        
    Tij = x(Nsq+N*i+j);
    Tijp1 = x(Nsq+N*i+j+1);
    Tijm1 = x(Nsq+N*i+j+1);
    Tip1j = x(Nsq+N*(i-1)+j);
    Tim1j = x(Nsq+N*(i-1)+j);
    Eijm1 = Eijp1+3.*h*sigma(j*h,i*h,Tij)*(4.-Eij);
    Eip1jm1 = Eip1jp1+3.*h*sigma(j*h,(i+1)*h,Tip1j)*(4.-Eip1j);
    Eim1jm1 = Eim1jp1+3.*h*sigma(j*h,(i-1)*h,Tim1j)*(4.-Eim1j);
    Eijp05 = (Eij+Eijp1)/2.;
    Eijm05 = (Eij+Eijm1)/2.;
    Eip05j = (Eij+Eip1j)/2.;
    Eim05j = (Eij+Eim1j)/2.;
    Tijp05 = (Tij+Tijp1)/2.;
    Tijm05 = (Tij+Tijm1)/2.;
    Tip05j = (Tij+Tip1j)/2.;
    Tim05j = (Tij+Tim1j)/2.;
    gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
    gradEijp05 = sqrt(gradEijp05);
    gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
    gradEijm05 = sqrt(gradEijm05);
    gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
    gradEip05j = sqrt(gradEip05j);
    gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
    gradEim05j = sqrt(gradEim05j);

    // computing diffusion
    fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                        -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                        +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                        -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
    fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                            -D2(Tijm05)*(Tij  -Tijm1)
                            +D2(Tip05j)*(Tip1j-Tij)
                            -D2(Tim05j)*(Tij  -Tim1j));

    //adding radiation
    radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
    fx(i*N+j) += radiation;
    fx(Nsq+i*N+j) -= radiation;
    
    
    //lower left corner
    i=0;
    j=0;
    //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=1 for E
    //yielding that value for Eijm1
    //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
    Eij = x(N*i+j);
    Eijp1 = x(N*i+j+1);        
    Eip1j = x(N*(i+1)+j);
    Eim1j = x(N*(i+1)+j);
    Eip1jp1 = x(N*(i+1)+j+1);
    Eim1jp1 = x(N*(i+1)+j+1);        
    Tij = x(Nsq+N*i+j);
    Tijp1 = x(Nsq+N*i+j+1);
    Tijm1 = x(Nsq+N*i+j+1);
    Tip1j = x(Nsq+N*(i+1)+j);
    Tim1j = x(Nsq+N*(i+1)+j);
    Eijm1 = Eijp1+3.*h*sigma(j*h,i*h,Tij)*(4.-Eij);
    Eip1jm1 = Eip1jp1+3.*h*sigma(j*h,(i+1)*h,Tip1j)*(4.-Eip1j);
    Eim1jm1 = Eim1jp1+3.*h*sigma(j*h,(i-1)*h,Tim1j)*(4.-Eim1j);
    Eijp05 = (Eij+Eijp1)/2.;
    Eijm05 = (Eij+Eijm1)/2.;
    Eip05j = (Eij+Eip1j)/2.;
    Eim05j = (Eij+Eim1j)/2.;
    Tijp05 = (Tij+Tijp1)/2.;
    Tijm05 = (Tij+Tijm1)/2.;
    Tip05j = (Tij+Tip1j)/2.;
    Tim05j = (Tij+Tim1j)/2.;
    gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
    gradEijp05 = sqrt(gradEijp05);
    gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
    gradEijm05 = sqrt(gradEijm05);
    gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
    gradEip05j = sqrt(gradEip05j);
    gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
    gradEim05j = sqrt(gradEim05j);

    // computing diffusion
    fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                        -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                        +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                        -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
    fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                            -D2(Tijm05)*(Tij  -Tijm1)
                            +D2(Tip05j)*(Tip1j-Tij)
                            -D2(Tim05j)*(Tij  -Tim1j));

    //adding radiation
    radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
    fx(i*N+j) += radiation;
    fx(Nsq+i*N+j) -= radiation;
    
    
    //upper right corner
    //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=0 for E
    //yielding that value for Eijp1
    i=N-1;
    j=N-1;
    //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
    Eij = x(N*i+j);        
    Eijm1 = x(N*i+j-1);
    Eip1j = x(N*(i-1)+j);
    Eim1j = x(N*(i-1)+j);
    Eip1jm1 = x(N*(i-1)+j-1);
    Eim1jm1 = x(N*(i-1)+j-1);
    Tij = x(Nsq+N*i+j);
    Tijp1 = x(Nsq+N*i+j-1);
    Tijm1 = x(Nsq+N*i+j-1);
    Tip1j = x(Nsq+N*(i-1)+j);
    Tim1j = x(Nsq+N*(i-1)+j);
    Eijp1 = Eijm1-3.*h*sigma(j*h,i*h,Tij)*Eij;
    Eip1jp1 = Eip1jm1-3.*h*sigma(j*h,(i+1)*h,Tip1j)*Eip1j;
    Eim1jp1 = Eim1jm1-3.*h*sigma(j*h,(i-1)*h,Tim1j)*Eim1j;
    Eijp05 = (Eij+Eijp1)/2.;
    Eijm05 = (Eij+Eijm1)/2.;
    Eip05j = (Eij+Eip1j)/2.;
    Eim05j = (Eij+Eim1j)/2.;
    Tijp05 = (Tij+Tijp1)/2.;
    Tijm05 = (Tij+Tijm1)/2.;
    Tip05j = (Tij+Tip1j)/2.;
    Tim05j = (Tij+Tim1j)/2.;
    gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
    gradEijp05 = sqrt(gradEijp05);
    gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
    gradEijm05 = sqrt(gradEijm05);
    gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
    gradEip05j = sqrt(gradEip05j);
    gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
    gradEim05j = sqrt(gradEim05j);

    // computing diffusion
    fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                        -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                        +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                        -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
    fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                            -D2(Tijm05)*(Tij  -Tijm1)
                            +D2(Tip05j)*(Tip1j-Tij)
                            -D2(Tim05j)*(Tij  -Tim1j));

    //adding radiation
    radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
    fx(i*N+j) += radiation;
    fx(Nsq+i*N+j) -= radiation;
    
    //lower right corner
    //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=0 for E
    //yielding that value for Eijp1
    i=0;
    j=N-1;
    //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
    Eij = x(N*i+j);        
    Eijm1 = x(N*i+j-1);
    Eip1j = x(N*(i+1)+j);
    Eim1j = x(N*(i+1)+j);
    Eip1jm1 = x(N*(i+1)+j-1);
    Eim1jm1 = x(N*(i+1)+j-1);
    Tij = x(Nsq+N*i+j);
    Tijp1 = x(Nsq+N*i+j-1);
    Tijm1 = x(Nsq+N*i+j-1);
    Tip1j = x(Nsq+N*(i+1)+j);
    Tim1j = x(Nsq+N*(i+1)+j);
    Eijp1 = Eijm1-3.*h*sigma(j*h,i*h,Tij)*Eij;
    Eip1jp1 = Eip1jm1-3.*h*sigma(j*h,(i+1)*h,Tip1j)*Eip1j;
    Eim1jp1 = Eim1jm1-3.*h*sigma(j*h,(i-1)*h,Tim1j)*Eim1j;
    Eijp05 = (Eij+Eijp1)/2.;
    Eijm05 = (Eij+Eijm1)/2.;
    Eip05j = (Eij+Eip1j)/2.;
    Eim05j = (Eij+Eim1j)/2.;
    Tijp05 = (Tij+Tijp1)/2.;
    Tijm05 = (Tij+Tijm1)/2.;
    Tip05j = (Tij+Tip1j)/2.;
    Tim05j = (Tij+Tim1j)/2.;
    gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
    gradEijp05 = sqrt(gradEijp05);
    gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
    gradEijm05 = sqrt(gradEijm05);
    gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
    gradEip05j = sqrt(gradEip05j);
    gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
    gradEim05j = sqrt(gradEim05j);

    // computing diffusion
    fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
                        -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
                        +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
                        -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
    fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
                            -D2(Tijm05)*(Tij  -Tijm1)
                            +D2(Tip05j)*(Tip1j-Tij)
                            -D2(Tim05j)*(Tij  -Tim1j));

    //adding radiation
    radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
    fx(i*N+j) += radiation;
    fx(Nsq+i*N+j) -= radiation;
}

Real RadiationDiffusion::D1(Real x1, Real x2, Real E, Real gradE, Real T)
{
    return 1./(3.*sigma(x1,x2,T)+gradE/E);
}

Real RadiationDiffusion::D2(Real T)
{
    return k*pow(T,2.5);
}

Real RadiationDiffusion::sigma(Real x1, Real x2, Real T)
{
    Real Z = (abs(x1-0.5)<1./6. && abs(x2-0.5)<1./6.) ? 10.:1.;
       
    return Z*Z*Z/T/T/T;
}

IntegroDifferentialEquation::IntegroDifferentialEquation()
{
    problem_name = "IntegroDifferentialEquation";
    
    tend=1.0;
    
    cte_rho = false;
    know_rho = false;
    analytical_df=true;
    dense_Jacobian=true;
    
    N = 100; 
    neqn = N;
    
    h = 1./N;
    
    sigma = 0.01;
}

void IntegroDifferentialEquation::set_initial_value(Vector& y0)
{
    const Real pi=4.*atan(1.);
    
    for(unsigned int j=0;j<N;j++)
        y0(j)= pow(cos((j+1.)*h*pi/2.),2);
}

void IntegroDifferentialEquation::f(Real t, Vector& x, Vector& fx)
{
    int i;
    Real xi;
    
    fx(0) = N*N*(x(1)-2.*x(0)+1.-sqrt(t)/2.); // diff with dirichlet 1-sqrt(t)/2
    for(i=1;i<N-1;i++)
    {
        fx(i) = N*N*(x(i-1)-2.*x(i)+x(i+1));
        
        xi = (i+1.)*h;
        fx(i) -= sigma*h/2.*( pow(1.-sqrt(t)/2.,4)/pow(1.+xi,2) + pow(x(N-1),4)/pow(2.-xi,2) );
        for(int j=1;j<N;j++)
            fx(i) -= sigma*h*pow(x(j-1),4)/pow(1.+abs(xi-j*h),2);
    }
    fx(N-1) = N*N*2.*(x(N-2)-x(N-1)); // diff with neumann


    i=0;
    xi = (i+1.)*h;
    fx(i) -= sigma*h/2.*( pow(1.-sqrt(t)/2.,4)/pow(1.+xi,2) + pow(x(N-1),4)/pow(2.-xi,2) );
    for(int j=1;j<N;j++)
        fx(i) -= sigma*h*pow(x(j-1),4)/pow(1.+abs(xi-j*h),2);
    
    i=N-1;
    xi = (i+1.)*h;
    fx(i) -= sigma*h/2.*( pow(1.-sqrt(t)/2.,4)/pow(1.+xi,2) + pow(x(N-1),4)/pow(2.-xi,2) );
    for(int j=1;j<N;j++)
        fx(i) -= sigma*h*pow(x(j-1),4)/pow(1.+abs(xi-j*h),2);
    
}

void IntegroDifferentialEquation::AN_df(Real t, Vector& x, Matrix& dfx)
{
    for(unsigned int i=0;i<N;i++)
    {
        for(unsigned int j=0;j<N;j++)
            dfx(i,j) = 4.*h*pow(x(j),3)/pow( 1+abs((i+1)*h-(j+1)*h) ,2);
        dfx(i,N-1) /= 2.;
    }
    
    dfx(0,0) -= 2.*N*N;
    dfx(0,1) += 1.*N*N;
    for(unsigned int i=1;i<N-1;i++)
    {
        dfx(i,i) -= 2.*N*N;
        dfx(i,i+1) += 1.*N*N;
        dfx(i,i-1) += 1.*N*N;
    }
    dfx(N-1,N-1) -= 2.*N*N;
    dfx(N-1,N-2) += 2.*N*N;
}



MonoDomain::MonoDomain()
{
    problem_name = "MonoDomain";
    
//    im = make_unique<uBidomain::HodgkinHuxley>();
    im = make_unique<uBidomain::MyoKit::MyoKitIonicModel>();
    
//    tend=800.; // [ms]
    tend = im->get_tend();        
    
    string unit = "mm";
    
    if(unit==string("cm"))
    {
        Lx = 2.; // cm
        Ly = 0.7; // cm
        Real h = 0.025; // cm //aimed mesh size, will be adapted a bit

        Nx = ceil(1+Lx/h);
        Ny = ceil(1+Ly/h);

        hx = Lx/(Nx-1.); //actual mesh size
        hy = Ly/(Ny-1.);

        init_mesh_data();    

        N_gating_var = im->get_n_gating_vars();
        neqn = (1+N_gating_var)*Nx*Ny;        

        cte_rho = false;
        know_rho = false;
        analytical_df=false;
        dense_Jacobian = false;

        diff_vect.resize(Nx*Ny);

        Real si_l = 1.7; //mS/cm
        Real se_l = 6.2; //mS/cm
        Real si_t = 0.19;
        Real se_t = 2.4;
        D_l = si_l*se_l/(si_l+se_l);
        D_t = si_t*se_t/(si_t+se_t);

        //using effective De for a cable equation
    //    Real k = 0.11;
    //    Real L = 125*1e-4;
    //    D = 1./(1./si+1./L/k);

        Cm = 1.; // uF/cm^2
        beta = 1400.; // cm^-1
        
        // im->Iion is in uA/cm^2, hence if we use cm there is no need to scale it
        scale_Iion = 1.;

        // define the stimulus region
        stim_vect = Vector::Zero(Nx*Ny);
        Real xmin = 0.;
        Real xmax = 0.16;
        Real ymin = 0.;
        Real ymax = 0.16;
        Real x,y;
        for(unsigned i=0;i<Ny;i++)
        for(unsigned j=0;j<Nx;j++)
        {
            x = j*hx;
            y = i*hy;
            if( x<=xmax && x>=xmin && y<=ymax && y>=ymin)
                stim_vect(i*Nx+j)=1.;
        }
    }
    else
    {
        Lx = 20.; // mm
        Ly = 7; // mm
        Real h = 0.25; // mm //aimed mesh size, will be adapted a bit

        Nx = ceil(1+Lx/h);
        Ny = ceil(1+Ly/h);

        hx = Lx/(Nx-1.); //actual mesh size
        hy = Ly/(Ny-1.);

        init_mesh_data();    

        N_gating_var = im->get_n_gating_vars();
        neqn = (1+N_gating_var)*Nx*Ny;        

        cte_rho = false;
        know_rho = false;
        analytical_df=false;
        dense_Jacobian = false;

        diff_vect.resize(Nx*Ny);

        Real si_l = 0.17; //mS/mm
        Real se_l = 0.62; //mS/mm
        Real si_t = 0.019;
        Real se_t = 0.24;
        D_l = si_l*se_l/(si_l+se_l);
        D_t = si_t*se_t/(si_t+se_t);

        Cm = 0.01; // uF/mm^2
        beta = 140.; // mm^-1
        
        // im->Iion is in uA/cm^2, hence if we use mm we scale it to
        // get a value in uA/mm^2
        scale_Iion = 1e-2;

        // define the stimulus region
        stim_vect = Vector::Zero(Nx*Ny);
        Real xmin = 0.;
        Real xmax = 1.6;
        Real ymin = 0.;
        Real ymax = 1.6;
        Real x,y;
        for(unsigned i=0;i<Ny;i++)
        for(unsigned j=0;j<Nx;j++)
        {
            x = j*hx;
            y = i*hy;
            if( x<=xmax && x>=xmin && y<=ymax && y>=ymin)
                stim_vect(i*Nx+j)=1.;
        }
    }
}

void MonoDomain::set_initial_value(Vector& y0)
{    
    y0.resize(neqn);
    im->get_y0(y0);    
}

void MonoDomain::init_mesh_data()
{
    points.resize(2*Nx*Ny);
    elements.resize(4*(Nx-1)*(Ny-1));
    
    for(unsigned int i=0;i<Ny;i++)
    for(unsigned int j=0;j<Nx;j++)
    {
        points[2*(i*Nx+j)] = j*hx;
        points[2*(i*Nx+j)+1] = i*hy;
    }
    
    for(unsigned int i=0;i<Ny-1;i++)
    for(unsigned int j=0;j<Nx-1;j++)
    {
        elements[4*(i*(Nx-1)+j)] = i*Nx+j;
        elements[4*(i*(Nx-1)+j)+1] = i*Nx+j+1;
        elements[4*(i*(Nx-1)+j)+2] = (i+1)*Nx+j+1;
        elements[4*(i*(Nx-1)+j)+3] = (i+1)*Nx+j;
    }
}

void MonoDomain::f(Real t, Vector& y, Vector& fy)
{   
    const auto& Vn = y.head(Nx*Ny);
    const auto& zn = y.tail(N_gating_var*Nx*Ny);    
    
    fy = im->f(t,y,false);    
    
    if(scale_Iion!=1.)
        fy.head(Nx*Ny) *= scale_Iion;
    fy.head(Nx*Ny) += diff(Vn)/beta+im->Iapp(t)*stim_vect*scale_Iion;
    
    if(Cm!=1.0)
        fy.head(Nx*Ny) /= Cm;
}

Vector MonoDomain::diff(const Vector& y)
{   
    Real facx = D_l/hx/hx;
    Real facy = D_t/hy/hy;
    
    diff_vect(0)     = 2.*facx*(y(1)-y(0))           //Lower-Left corner
                      +2.*facy*(y(Nx)-y(0));
    diff_vect(Nx-1)   = 2.*facx*(y(Nx-2)-y(Nx-1))          //Lower-Right corner
                      +2.*facy*(y(2*Nx-1)-y(Nx-1));
    diff_vect(Nx*Ny-Nx) = 2.*facx*(y(Nx*Ny-Nx+1)-y(Nx*Ny-Nx))       //Upper-Left corner
                      +2.*facy*(y(Nx*Ny-2*Nx)-y(Nx*Ny-Nx));
    diff_vect(Nx*Ny-1) = 2.*facx*(y(Nx*Ny-2)-y(Nx*Ny-1))            //Upper-Right corner
                      +2.*facy*(y(Nx*Ny-Nx-1)-y(Nx*Ny-1));
    
    unsigned i,j;
    
    #pragma omp parallel for schedule(dynamic) private(j) shared(Nx,Ny,diff_vect,facx,facy,y)
    for(j=1;j<Nx-1;j++)
    {
        diff_vect(j) = facx*(y(j-1)-2.*y(j)+y(j+1))           // Lower boundary
                       +2.*facy*((y(j+Nx)-y(j)));
        diff_vect(Nx*Ny-Nx+j) = facx*(y(Nx*Ny-Nx+j-1)-2.*y(Nx*Ny-Nx+j)+y(Nx*Ny-Nx+j+1))   // Upper boundary
                               +2.*facy*((y(Nx*Ny-2*Nx+j)-y(Nx*Ny-Nx+j)));
    }
    #pragma omp parallel for schedule(dynamic) private(i) shared(Nx,Ny,diff_vect,facx,facy,y)
    for(i=1;i<Ny-1;i++)
    {
        diff_vect(i*Nx) = 2.*facx*(y(i*Nx+1)-y(i*Nx))                 // Left boundary
                          +facy*(y(i*Nx-Nx)-2.*y(i*Nx)+y(i*Nx+Nx));
        diff_vect(i*Nx+Nx-1) = 2.*facx*(y(i*Nx+Nx-2)-y(i*Nx+Nx-1))          // Right boundary
                             +facy*(y(i*Nx+2*Nx-1)-2.*y(i*Nx+Nx-1)+y(i*Nx-1));
    }
    
    #pragma omp parallel for schedule(dynamic) private(i,j) shared(Nx,Ny,diff_vect,facx,facy,y)
    for(i=1;i<Ny-1;i++)
    for(j=1;j<Nx-1;j++)
        diff_vect(i*Nx+j) = facx*(y(i*Nx+j+1)-2.*y(i*Nx+j)+y(i*Nx+j-1))       // Internal points
                           +facy*(y(i*Nx+j+Nx)-2.*y(i*Nx+j)+y(i*Nx+j-Nx));
   
    return diff_vect;
}


void MonoDomain::write_solution(const int nout, const string solname, 
                                const Real t, const Vector& y)
{
    const Vector& u = y.segment(0,Nx*Ny);
    vector<Real> u_stdvec(u.data(), u.data() + u.size());
    
    // write a VTK solution
    const int dim = 2;
    const int cell_size = 4;
   
    std::string extension = ".vtu";
    std::string filename_u = solname+"_u_"+to_string(nout)+extension;

    leanvtk::VTUWriter writer;
    
    writer.add_scalar_field("u", u_stdvec);
    
    for(unsigned int k=1;k<=N_gating_var;k++)
    {
        const Vector& g = y.segment(Nx*Ny*k,Nx*Ny);
        vector<Real> g_stdvec(g.data(), g.data() + g.size());
        writer.add_scalar_field(string("g_")+to_string(k), g_stdvec);
    }
    writer.write_surface_mesh(filename_u, dim, cell_size, points, elements);
}