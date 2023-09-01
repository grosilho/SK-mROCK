#include "Ode.h"

Ode::Ode()
{
    cte_rho=false;
    know_rho=false;
    analytical_df=false;
    autodiff_df = false;
    dense_Jacobian=true;
    
    multirate=false;
}

int Ode::neqn;

Ode::~Ode()
{
}

void Ode::rho(Real t, Vector& y, Real& eigmax)
{
    cout<<"ERROR: using a non implemented rho function!\nUse -intrho 1."<<endl;
}

bool Ode::estimation_rho_known()
{
    return know_rho;
}

bool Ode::analytical_df_known()
{
    return analytical_df;
}

Real Ode::get_tend()
{
    return tend;
}

int Ode::get_system_size()
{
    return neqn;
}

bool Ode::is_rho_constant()
{
    return cte_rho;
}

string Ode::get_problem_name() 
{
    return problem_name;
}

bool Ode::is_multirate()
{
    return multirate;
}

bool Ode::has_dense_Jacobian()
{
    return dense_Jacobian;
}

bool Ode::autodiff_implemented()
{
    return autodiff_df;
}

void Ode::df(Real t, Vector& x, Matrix& dfx)
{
    if(analytical_df)
        this->AN_df(t,x,dfx);
    else if(autodiff_df)
        this->AD_df(t,x,dfx);
    else//finite difference
        this->FD_df(t,x,dfx);
}

void Ode::df(Real t, Vector& x, SpMatrix& dfx)
{
    if(analytical_df)
        this->AN_df(t,x,dfx);
    else if(autodiff_df)
        this->AD_df(t,x,dfx);
    else//finite difference
        this->FD_df(t,x,dfx);
}

void Ode::FD_df(Real t, Vector& x, Matrix& dfx)
{
    dfx.resize(neqn,neqn);
    Real epsilon = 1e-10;
    static Vector fx(neqn), fxtmp(neqn), xtmp(neqn);
    f(t,x,fx);
    xtmp = x;
    
    for(unsigned int i=0;i<neqn;i++)
    {
        xtmp(i) = xtmp(i) + epsilon*(1.+abs(xtmp(i)));
        f(t,xtmp,fxtmp);
        dfx.block(0,i,neqn,1) = (fxtmp-fx)/(epsilon*(1.+abs(xtmp(i))));
        xtmp(i)=x(i);
    }
}

void Ode::FD_df(Real t, Vector& x, SpMatrix& dfx)
{
    dfx.resize(neqn,neqn);
    dfx.reserve(Eigen::VectorXi::Constant(neqn,8));
    
    Real epsilon = 1e-10;
    static Vector fx(neqn), fxtmp(neqn), xtmp(neqn);
    f(t,x,fx);
    xtmp = x;
    Real tol_nz = 1e-10;
    
    for(unsigned int i=0;i<neqn;i++)
    {
        xtmp(i) = xtmp(i) + epsilon*(1.+abs(xtmp(i)));
        f(t,xtmp,fxtmp);
        fxtmp = (fxtmp-fx)/(epsilon*(1.+abs(xtmp(i))));
        xtmp(i)=x(i);
        for(unsigned int j=0;j<neqn;j++)
            if(abs(fxtmp(j))>tol_nz*(1.+abs(xtmp(i))))
                dfx.insert(j,i)=fxtmp(j);
    }
    dfx.makeCompressed();  
}

void Ode::AN_df(Real t, Vector& x, Matrix& fx)
{
     cout<<"ERROR: using a non implemented AN_df function!"<<endl;
}

void Ode::AD_df(Real t, Vector& x, Matrix& fx)
{
     cout<<"ERROR: using a non implemented AD_df function!"<<endl;
}

void Ode::AN_df(Real t, Vector& x, SpMatrix& fx)
{
     cout<<"ERROR: using a non implemented AN_df function!"<<endl;
}

void Ode::AD_df(Real t, Vector& x, SpMatrix& fx)
{
     cout<<"ERROR: using a non implemented AD_df function!"<<endl;
}


void Ode::write_solution(const int nout, const string solname, 
                         const Real t, const Vector& y)
{
    
}


// FAST SLOW ODE ---------------------------------------------------------------

MultirateOde::MultirateOde()
{
    multirate=true;
}

MultirateOde::~MultirateOde()
{
    
}

void MultirateOde::rho(Real t, Vector& y, 
                       Real& eigmax_F, Real& eigmax_S)
{
    cout<<"ERROR: using a non implemented rho function!"<<endl;
}

void MultirateOde::dfF(Real t, Vector& x, Matrix& dfx)
{
    dfx.resize(neqn,neqn);
    Real epsilon = 1e-10;
    static Vector fx(neqn), fxtmp(neqn), xtmp(neqn);
    fF(t,x,fx);
    xtmp = x;
    
    for(unsigned int i=0;i<neqn;i++)
    {
        xtmp(i) = xtmp(i) + epsilon*(1.+abs(xtmp(i)));
        fF(t,xtmp,fxtmp);
        dfx.block(0,i,neqn,1) = (fxtmp-fx)/(epsilon*(1.+abs(xtmp(i))));
        xtmp(i)=x(i);
    }
}

void MultirateOde::dfF(Real t, Vector& x, SpMatrix& dfx)
{
    dfx.resize(neqn,neqn);
    dfx.reserve(Eigen::VectorXi::Constant(neqn,8));
    
    Real epsilon = 1e-10;
    static Vector fx(neqn), fxtmp(neqn), xtmp(neqn);
    fF(t,x,fx);
    xtmp = x;
    Real tol_nz = 1e-10;
    
    for(unsigned int i=0;i<neqn;i++)
    {
        xtmp(i) = xtmp(i) + epsilon*(1.+abs(xtmp(i)));
        fF(t,xtmp,fxtmp);
        fxtmp = (fxtmp-fx)/(epsilon*(1.+abs(xtmp(i))));
        xtmp(i)=x(i);
        for(unsigned int j=0;j<neqn;j++)
            if(abs(fxtmp(j))>tol_nz*(1.+abs(xtmp(i))))
                dfx.insert(j,i)=fxtmp(j);
    }
    dfx.makeCompressed();  
}