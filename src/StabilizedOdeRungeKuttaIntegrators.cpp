#include "StabilizedOdeRungeKuttaIntegrators.h"
#include "ChebyshevMethods.h"
#include "LegendreMethods.h"

IERKC::IERKC(Parameters* param_, MultirateOde* mode_)
:MultirateOdeRungeKuttaIntegrator(param_,mode_)
{
    Astable = false;
    err_order = 1;
    
    damping = 0.05;
    beta = 2.-4.*damping/3.;
       
    for(int i=0;i<4;i++)
        integr_add[i]=new Vector(mode->get_system_size());
    
    if(ode->has_dense_Jacobian())
    {
        J.resize(ode->get_system_size(),ode->get_system_size());
        I = Eigen::MatrixXd::Identity(ode->get_system_size(),ode->get_system_size());
    }
    else//use sparse matrices
    {
        spJ.resize(ode->get_system_size(),ode->get_system_size());
        spI.resize(ode->get_system_size(),ode->get_system_size());
        spI.setIdentity();
    }
    
    reinit_integrator();
}

IERKC::~IERKC()
{
    for(int i=0;i<4;i++)
        delete integr_add[i];
}

void IERKC::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by RKC. It follows the formula
     * \f$\rho \Delta t \leq \beta s^2 \f$ with some cautions.
     */
    
    unsigned int stages_limit = 1e6;
    unsigned int safe_add = 1;
    
    if(h*eigmax_S<1.5) //for s=1 we have beta=2 and not 2-4/3*damping, so maybe one stage is enough
        s=1;
    else if(h*eigmax_S<2.) // if close to 2 then add some stages for safety
        s=1+safe_add;
    else // if EE not stable then use RKC formula
    {
        s = ceil(sqrt(h*eigmax_S/beta));
        s += safe_add;
    }

    //we try to avoid more than max_s stages for roundoff errors reasons
    if(s>stages_limit)
    {
        //if more than max_s stages are needed we resize dt
        s=stages_limit; 
        h=0.95*beta*s*s/eigmax_S;     
        last=false;
        err_control.update_hn(h);
    }
    else if(dt_adaptivity)
        s = max(s,2); // at least two stages, needed to estimate the error 
    
//    write_parameters(h);
    
    s_max= max(s_max,s);
    s_avg += s;    
}

void IERKC::write_parameters(Real dt)
{
    static int nout=1;
    ofstream outfile;
            
    if(nout==1)
    {
        outfile.open(param->output_path+string("_params.csv"), ofstream::out);
        outfile<<"s, dt, rhoF, rhoS"<<endl;
    }
    else
        outfile.open(param->output_path+string("_params.csv"), ofstream::out | ofstream::app);
    
    outfile<<setprecision(16)<<s<<", "<<dt<<", "<<eigmax_F<<", "<<eigmax_S<<endl;      
     
    outfile.close();
    
    nout++;
}

void IERKC::step(const Real t, const Real& h)
{
    Vector*& fn= integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> d;
    vector<Real> nu;
    vector<Real> mu;
    vector<Real> kappa;
    
    ChebyshevMethods::CoefficientsRKC1(mu, nu, kappa, c, d, s, damping);
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    mode->fS(t,*Kjm1,*fn);     
    *Kj = mu[0]*h*(*fn);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        mode->fS(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);                                   
    }   
    
    static unsigned nstep=0;
    
    Vector*& dk= integr[0];
    Vector*& hfyk= integr[1];
    Vector*& ddk= integr[2];
    Vector*& tmp1= integr[3];
    
            
    if(ode->has_dense_Jacobian())
    {
        if(nstep==0)
        {
            nstep++;
            mode->dfF(t+h,*Kj,J);
            J = I - J*h;
            dir_solver.compute(J);
        }
        *ynpu = dir_solver.solve(*Kj);
    }
    else
    {
        if(nstep==0)
        {
            nstep++;
            mode->dfF(t+h,*Kj,spJ);
            spJ = spI - spJ*h;
            solver.compute(spJ);
        }
        *ynpu = solver.solve(*Kj);
        lin_solv_iter = solver.iterations();
    }
    
    n_fS_eval += s;
    n_fF_eval += 1;
}  

void IERKC::disp_step_info(Real& t, Real& h, bool accepted)
{
    std::string delta = u8"\u0394";
    string rho =u8"\u03C1";
    
    cout << scientific;
    
    cout<<"Step t = "<<setw(6)<<setprecision(4)<<t<<", "<<delta<<"t = "<<setw(8)<<setprecision(6)<<h
    <<", s = "<<setw(3)<<s<<", m = "<<setw(3)<<m
    <<", "<<rho<<"F = "<<setw(3)<<eigmax_F<<", "<<rho<<"S = "<<setw(3)<<eigmax_S
    <<" and |y_n| = "<<setw(7)<<setprecision(4)<<ynpu->lpNorm<Eigen::Infinity>()<<". ";
    if(accepted)
        cout<<" Accepted ";
    else
        cout<<" Rejected ";
    err_control.disp_info();
    cout<<endl;
}

// ----------------------------------------------------------------------------

mRKC::mRKC(Parameters* param_, MultirateOde* mode_)
:MultirateOdeRungeKuttaIntegrator(param_,mode_)
{
    Astable = false;
    err_order = 1;
    
    damping = 0.05;
    beta = 2.-4.*damping/3.;
       
    for(int i=0;i<4;i++)
        integr_add[i]=new Vector(mode->get_system_size());
    
    reinit_integrator();
}

mRKC::~mRKC()
{
    for(int i=0;i<4;i++)
        delete integr_add[i];
}

void mRKC::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by RKC. It follows the formula
     * \f$\rho \Delta t \leq \beta s^2 \f$ with some cautions.
     */
    
    unsigned int stages_limit = 1e6;
    unsigned int safe_add = 1;
    

    if(h*eigmax_S<1.5) //for s=1 we have beta=2 and not 2-4/3*damping, so maybe one stage is enough
        s=1;
    else if(h*eigmax_S<2.) // if close to 2 then add some stages for safety
        s=1+safe_add;
    else // if EE not stable then use RKC formula
    {
        s = ceil(sqrt(h*eigmax_S/beta));
        s += safe_add;
    }

    //we try to avoid more than max_s stages for roundoff errors reasons
    if(s>stages_limit)
    {
        //if more than max_s stages are needed we resize dt
        s=stages_limit; 
        h=0.95*beta*s*s/eigmax_S;     
        last=false;
        err_control.update_hn(h);
    }
    else if(dt_adaptivity)
        s = max(s,2); // at least two stages, needed to estimate the error

    // General ODE
//    m = ceil(sqrt(6.*eigmax_F*h/beta/beta/s/s+1.))+safe_add;
//    eta = 6.*m*m*h/(m*m-1.)/beta/s/s;
    
    // Parabolic ODE in refined mesh or ODE with scale separation
    eta = 2.*h/beta/s/s;
    m = ceil(sqrt(eta*eigmax_F/beta))+safe_add;    
    
//    write_parameters(h);
    
    s_max= max(s_max,s);
    m_max= max(m_max,m);
    s_avg += s;
    m_avg += m;
    
}

void mRKC::write_parameters(Real dt)
{
    static int nout=1;
    ofstream outfile;
            
    if(nout==1)
    {
        outfile.open(param->output_path+string("_params.csv"), ofstream::out);
        outfile<<"s, m, dt, eta, rhoF, rhoS"<<endl;
    }
    else
        outfile.open(param->output_path+string("_params.csv"), ofstream::out | ofstream::app);
    
    outfile<<setprecision(16)<<s<<", "<<m<<", "<<dt<<", "<<eta<<", "<<eigmax_F<<", "<<eigmax_S<<endl;      
     
    outfile.close();
    
    nout++;
}

void mRKC::step(const Real t, const Real& h)
{
    Vector*& fn= integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> d;
    vector<Real> nu;
    vector<Real> mu;
    vector<Real> kappa;
    
    ChebyshevMethods::CoefficientsRKC1(mu, nu, kappa, c, d, s, damping);
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    f_eta(t,*Kjm1,*fn); 
    *Kj = mu[0]*h*(*fn);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        f_eta(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);                                   
    }   
    
    swap_ptr=ynpu;
    ynpu=Kj;
    Kj=swap_ptr;

    //compute errors
    if(s>=2)
    {
        Real mult, a1, a2, a3;    
        if(s>2)
        {
            mult = (d[s-1]-1.)/((c[s-3]-1.)*(d[s-2]-d[s-1])-(c[s-2]-1)*(d[s-3]-d[s-1]));
            a1 = (1.-c[s-2])*mult;
            a2 = (c[s-3]-1.)*mult;
            a3 = (c[s-2]-c[s-3])*mult;
        }
        else
        {        
            mult = (d[s-1]-1.)/d[s-1]/c[s-2];
            a1 = (1.-c[s-2])*mult;
            a2 = -mult;
            a3 = c[s-2]*mult;
        }
        *Kjm2 *= a1;
        *Kjm2 += a2*(*Kjm1);
        *Kjm2 += a3*(*ynpu);

        for(int i=0;i<mode->get_system_size();i++)
            (*Kjm2)(i) = (*Kjm2)(i)/(a_tol+r_tol*max(abs((*ynpu)(i)),abs((*yn)(i))));

        err_control.get_errD() =Kjm2->norm()/sqrt(mode->get_system_size());
    }
    else
        err_control.get_errD() = 0.;
    
    n_fS_eval += s;
    n_fF_eval += s*m;
        
}  

void mRKC::f_eta(Real t, Vector& x, Vector& fx)
{
    Vector*& Ljm1= integr_add[0];
    Vector*& Ljm2= integr_add[1];
    Vector*& Lj= integr_add[2];
    Vector*& R0= integr_add[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> d;
    vector<Real> nu;
    vector<Real> mu;
    vector<Real> kappa;
    
    ChebyshevMethods::CoefficientsRKC1(mu, nu, kappa, c, d, m, damping);
    
    
    // Stage 0
    mode->fS(t,x,*R0);
    *Ljm1 = x;
    
    // Stage 1
    mode->fF(t,x,*Lj);
    *Lj += *R0;
    *Lj *= eta*mu[0];
    *Lj += *Ljm1;

    
    // Stages j=2,...,m
    for(unsigned int j=2;j<=m;j++)
    {
        swap_ptr=Ljm2;
        Ljm2=Ljm1;
        Ljm1=Lj;
        Lj=swap_ptr;               
                
        mode->fF(t+eta*c[j-2],*Ljm1,*Lj);   
        *Lj += *R0;
        *Lj *= eta*mu[j-1];
        *Lj += nu[j-1]*(*Ljm1);
        *Lj += kappa[j-1]*(*Ljm2);                                   
    }   
    
    fx = *Lj;
    fx -= x;
    fx *= (1./eta);        
}  

void mRKC::disp_step_info(Real& t, Real& h, bool accepted)
{
    std::string delta = u8"\u0394";
    string rho =u8"\u03C1";
    
    cout << scientific;
    
    cout<<"Step t = "<<setw(6)<<setprecision(4)<<t<<", "<<delta<<"t = "<<setw(8)<<setprecision(6)<<h
    <<", s = "<<setw(3)<<s<<", m = "<<setw(3)<<m<<", eta = "<<setw(3)<<eta
    <<", "<<rho<<"F = "<<setw(3)<<eigmax_F<<", "<<rho<<"S = "<<setw(3)<<eigmax_S
    <<" and |y_n| = "<<setw(7)<<setprecision(4)<<ynpu->lpNorm<Eigen::Infinity>()<<". ";
    if(accepted)
        cout<<" Accepted ";
    else
        cout<<" Rejected ";
    err_control.disp_info();
    cout<<endl;
}

// ----------------------------------------------------------------------------

RKC1::RKC1(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_,ode_)
{
    Astable = false;
    err_order = 1;
    
    damping = 0.1;
    beta = 2.-4.*damping/3.;

    reinit_integrator();
}

RKC1::~RKC1()
{

}

void RKC1::reinit_integrator()
{
    n_output_eigs=0;
    OdeRungeKuttaIntegrator::reinit_integrator();
}

void RKC1::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by RKC. It follows the formula
     * \f$\rho \Delta t \leq \beta s^2 \f$ with some cautions.
     */
    
    unsigned int stages_limit = 1e6;
    unsigned int safe_add = 2;

    if(h*eigmax<1.5)
        s = 1;
    else if(h*eigmax<2.) //for s=1 we have beta=2 and not 2-4/3*damping, so maybe one stage is enough
        s = 1+safe_add;
    else
    {
        if(damping<=0.1)
            s = ceil(sqrt(h*eigmax/beta));
        else
        {
            s = ceil(sqrt(h*eigmax/2.0));
            while(h*eigmax >= ChebyshevMethods::ls_RKC1(s,damping))
                ++s;
        }
        s += safe_add;
    }

    //we try to avoid more than max_s stages for roundoff errors reasons
    if(s>stages_limit)
    {
        //if more than max_s stages are needed we resize dt
        s=stages_limit; 
        if(damping<=0.1)
            h=0.9*beta*s*s/eigmax;     
        else
            h = 0.9*ChebyshevMethods::ls_RKC1(s,damping)/eigmax;
        last=false;
        err_control.update_hn(h);
    }
    else if(dt_adaptivity)
        s = max(s,2); // at least two stages, needed to estimate the error
            
    if(s>s_max)
        s_max=s;    
    s_avg += s;
    
//    write_eigenvalues();
}

void RKC1::write_eigenvalues()
{
    Matrix jac;
    ode->df(t,*yn,jac);
    Eigen::VectorXcd eigvals = jac.eigenvalues();
    
    
//    if(n_output_eigs==352)
//    {
//        ofstream outfile_mat;
//        outfile_mat.open(param->output_path+string("_mat.m"), ofstream::out);
//        outfile_mat<<setprecision(16);
//        outfile_mat<<"t="<<t<<";"<<endl;
//        outfile_mat<<"jac = zeros("<<jac.rows()<<","<<jac.cols()<<");"<<endl;
//        for(unsigned int j=0;j<jac.cols();j++)
//        {
//            outfile_mat<<"jac(:,"<<j+1<<")=[";
//            for(unsigned int i=0;i<jac.rows()-1;i++)
//              outfile_mat<<jac(i,j)<<"; ";
//            outfile_mat<<jac(jac.rows()-1,j)<<"];"<<endl;
//        }
//        outfile_mat.close();
//    }
    
    ofstream outfile;
    
    if(n_output_eigs==1)
        outfile.open(param->output_path+string("_eigs.m"), ofstream::out);
    else
        outfile.open(param->output_path+string("_eigs.m"), ofstream::out | ofstream::app);
    
    outfile<<setprecision(16);
    outfile<<"t("<<n_output_eigs<<")="<<t<<";"<<endl;
    outfile<<"eigvals("<<n_output_eigs<<",:)=[";
    for(auto e: eigvals.head(eigvals.size()-1))
        outfile<<e.real()<<"+"<<e.imag()<<"*i, ";
    auto e = eigvals.tail(1);
    outfile<<e.real()<<"+"<<e.imag()<<"*i];"<<endl;
    
    outfile.close();
    n_output_eigs++;
}

void RKC1::write_parameters(Real dt)
{
    static int nout=1;
    ofstream outfile;
            
    if(nout==1)
    {
        outfile.open(param->output_path+string("_params.csv"), ofstream::out);
        outfile<<"s, dt, rho"<<endl;
    }
    else
        outfile.open(param->output_path+string("_params.csv"), ofstream::out | ofstream::app);
    
    outfile<<setprecision(16)<<s<<", "<<dt<<", "<<eigmax<<endl;      
     
    outfile.close();
    
    nout++;
}

void RKC1::step(const Real t, const Real& h)
{
    Vector*& fn= integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> d;
    vector<Real> nu;
    vector<Real> mu;
    vector<Real> kappa;
    
    ChebyshevMethods::CoefficientsRKC1(mu, nu, kappa, c, d, s, damping);
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    ode->f(t,*Kjm1,*fn); 
    *Kj = mu[0]*h*(*fn);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        ode->f(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);                                   
    }   
    
    swap_ptr=ynpu;
    ynpu=Kj;
    Kj=swap_ptr;

    //compute errors
    Real mult, a1, a2, a3;    
    if(s>2)
    {
        mult = (d[s-1]-1.)/((c[s-3]-1.)*(d[s-2]-d[s-1])-(c[s-2]-1)*(d[s-3]-d[s-1]));
        a1 = (1.-c[s-2])*mult;
        a2 = (c[s-3]-1.)*mult;
        a3 = (c[s-2]-c[s-3])*mult;        
    }
    else// if(s==2)
    {        
        mult = (d[s-1]-1.)/d[s-1]/c[s-2];
        a1 = (1.-c[s-2])*mult;
        a2 = -mult;
        a3 = c[s-2]*mult;
    }
    *Kjm2 *= a1;
    *Kjm2 += a2*(*Kjm1);
    *Kjm2 += a3*(*ynpu);
     
    for(int i=0;i<ode->get_system_size();i++)
        (*Kjm2)(i) = (*Kjm2)(i)/(a_tol+r_tol*max(abs((*ynpu)(i)),abs((*yn)(i))));
    
    err_control.get_errD() =Kjm2->norm()/sqrt(ode->get_system_size());
            
    n_f_eval=n_f_eval+s;
}  

RKC2::RKC2(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_,ode_)
{
    Astable = false;
    err_order = 2;
    
    damping = 2./13.;
    beta = 2./3.*(1.-2.*damping/15.);
    
    reinit_integrator();
}

RKC2::~RKC2()
{

}

void RKC2::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by RKC. It follows the formula
     * \f$\rho \Delta t \leq 0.65 (s^2-1) \f$ with some cautions.
     */
    
    unsigned int stages_limit = 1e6;
    unsigned int safe_add = 2;
    
    if(damping<=0.1)
        s = ceil(sqrt(h*eigmax/beta+1.));
    else
    {
        s = ceil(sqrt(h*eigmax/(2./3.)+1.));
        while(h*eigmax >= ChebyshevMethods::ls_RKC2(s,damping))
                ++s;
    }    
    s += safe_add;
                
    //we try to avoid more than max_s stages for internal stability reasons
    if(s>stages_limit) 
    {
        //if more than max_s stages are needed we resize h       
        s = stages_limit;
        if(damping<=0.1)
            h = 0.9*beta*(s*s-1)/eigmax;
        else
            h = 0.9*ChebyshevMethods::ls_RKC2(s,damping)/eigmax;
        last=false;
        err_control.update_hn(h);
    }
    else//at least 2 stages
        s=max(s,2);
            
    if (s>s_max)
        s_max=s;
    s_avg += s;
}

void RKC2::step(const Real t, const Real& h)
{
    /**
     * Does a RKC step. Given \f$ y_n\f$ at \f$t\f$ it computes 
     * \f$y\f$ at \f$t+h\f$. Parabolic provides the right hand side function
     * and sets the Dirichlet boundary conditions.
     */
    
    Vector*& fn= integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> mu;
    vector<Real> nu;
    vector<Real> kappa;
    vector<Real> gamma;
    
    ChebyshevMethods::CoefficientsRKC2(mu, nu, kappa, gamma, c, s, damping);
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    ode->f(t,*Kjm1,*fn); 
    *Kj = mu[0]*h*(*fn);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        ode->f(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);     
        *Kj += gamma[j-1]*h*(*fn);
        *Kj += (1.-nu[j-1]-kappa[j-1])*(*yn);
    }   
    
    swap_ptr=ynpu;
    ynpu=Kj;
    Kj=swap_ptr;

    // computing local error
    ode->f(t+h,*ynpu,*Kjm1); 
    for(int i=0;i<ode->get_system_size();i++)
    {
        mu[0] = max(abs((*ynpu)(i)),abs((*yn)(i)));
        nu[0] = 0.8*((*yn)(i)-(*ynpu)(i))+0.4*h*((*fn)(i)+(*Kjm1)(i));
        (*Kjm2)(i) = nu[0]/(a_tol+mu[0]*r_tol);
    }
    err_control.get_errD() =Kjm2->norm()/sqrt(ode->get_system_size());
    
    n_f_eval=n_f_eval+s;
}  

// ----------------------------------------------------------------------------

RKL1::RKL1(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_,ode_)
{
    Astable = false;
    err_order = 1;
    
    damping = 0.0;
    // l_s = beta*s*(s+1);
    
    reinit_integrator();
}

RKL1::~RKL1()
{

}

void RKL1::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by RKC. It follows the formula
     * \f$\rho \Delta t \leq \beta s^2 \f$ with some cautions.
     */
    
    unsigned int stages_limit = 1e6;
    unsigned int safe_add = 2;

    if(h*eigmax<1.5)
        s = 1;
    else if(h*eigmax<2.) //for s=1 we have beta=2 and not depending on damping, so maybe one stage is enough
        s = 1+safe_add;
    else
    {
        s = ceil( (sqrt(1.+4.*h*eigmax)-1.)/2. );
        if(damping>0.)
            while(h*eigmax >= LegendreMethods::ls_RKL1(s,damping))
                ++s;
        s += safe_add;
    }
    
    //we try to avoid more than max_s stages for roundoff errors reasons
    if(s>stages_limit)
    {
        //if more than max_s stages are needed we resize dt
        s=stages_limit; 
        if(damping==0.)
            h=0.9*s*(s+1.)/eigmax;     
        else 
            h=0.9*LegendreMethods::ls_RKL1(s,damping)/eigmax;     
        last=false;
        err_control.update_hn(h);
    }
    else if(dt_adaptivity)
        s = max(s,2); // at least two stages, needed to estimate the error
            
    if(s>s_max)
        s_max=s;    
    s_avg += s;
}

void RKL1::step(const Real t, const Real& h)
{
    Vector*& fn= integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> d;
    vector<Real> nu;
    vector<Real> mu;
    vector<Real> kappa;
    
    LegendreMethods::CoefficientsRKL1(mu, nu, kappa, c, d, s, damping);
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    ode->f(t,*Kjm1,*fn); 
    *Kj = mu[0]*h*(*fn);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        ode->f(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);                                   
    }   
    
    swap_ptr=ynpu;
    ynpu=Kj;
    Kj=swap_ptr;

    //compute errors
    Real mult, a1, a2, a3;    
    if(s>2)
    {
        mult = (d[s-1]-1.)/((c[s-3]-1.)*(d[s-2]-d[s-1])-(c[s-2]-1)*(d[s-3]-d[s-1]));
        a1 = (1.-c[s-2])*mult;
        a2 = (c[s-3]-1.)*mult;
        a3 = (c[s-2]-c[s-3])*mult;        
    }
    else// if(s==2)
    {        
        mult = (d[s-1]-1.)/d[s-1]/c[s-2];
        a1 = (1.-c[s-2])*mult;
        a2 = -mult;
        a3 = c[s-2]*mult;
    }
    *Kjm2 *= a1;
    *Kjm2 += a2*(*Kjm1);
    *Kjm2 += a3*(*ynpu);
     
    for(int i=0;i<ode->get_system_size();i++)
        (*Kjm2)(i) = (*Kjm2)(i)/(a_tol+r_tol*max(abs((*ynpu)(i)),abs((*yn)(i))));
    
    err_control.get_errD() =Kjm2->norm()/sqrt(ode->get_system_size());
            
    n_f_eval=n_f_eval+s;
}  

RKL2::RKL2(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_,ode_)
{
    Astable = false;
    err_order = 2;
    
    damping = 0.;
    
    reinit_integrator();
}

RKL2::~RKL2()
{

}

void RKL2::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by RKC. It follows the formula
     * \f$\rho \Delta t \leq 0.65 (s^2-1) \f$ with some cautions.
     */
    
    unsigned int stages_limit = 1e6;
    unsigned int safe_add = 2;
    
    s = ceil(0.5*(sqrt(1.+8.*(1.+h*eigmax))-1.));
    if(damping>0.)
        while(h*eigmax >= LegendreMethods::ls_RKL2(s,damping))
            ++s;
    s+=safe_add;
    
    //we try to avoid more than max_s stages for internal stability reasons
    if(s>stages_limit) 
    {
        //if more than max_s stages are needed we resize h       
        s = stages_limit;
        if(damping==0.)
            h = 0.9*0.5*(s-1.)*(s+2.)/eigmax;
        else
            h = 0.9*LegendreMethods::ls_RKL2(s,damping)/eigmax;
        last=false;
        err_control.update_hn(h);
    }
    else//at least 2 stages
        s=max(s,2);
            
    if (s>s_max)
        s_max=s;
    s_avg += s;
}

void RKL2::step(const Real t, const Real& h)
{
    /**
     * Does a RKC step. Given \f$ y_n\f$ at \f$t\f$ it computes 
     * \f$y\f$ at \f$t+h\f$. Parabolic provides the right hand side function
     * and sets the Dirichlet boundary conditions.
     */
    
    Vector*& fn= integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> mu;
    vector<Real> nu;
    vector<Real> kappa;
    vector<Real> gamma;
    
    LegendreMethods::CoefficientsRKL2(mu, nu, kappa, gamma, c, s, damping);
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    ode->f(t,*Kjm1,*fn); 
    *Kj = mu[0]*h*(*fn);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        ode->f(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);     
        *Kj += gamma[j-1]*h*(*fn);
        *Kj += (1.-nu[j-1]-kappa[j-1])*(*yn);
    }   
    
    swap_ptr=ynpu;
    ynpu=Kj;
    Kj=swap_ptr;

    // computing local error
    ode->f(t+h,*ynpu,*Kjm1); 
    for(int i=0;i<ode->get_system_size();i++)
    {
        mu[0] = max(abs((*ynpu)(i)),abs((*yn)(i)));
        nu[0] = 0.8*((*yn)(i)-(*ynpu)(i))+0.4*h*((*fn)(i)+(*Kjm1)(i));
        (*Kjm2)(i) = nu[0]/(a_tol+mu[0]*r_tol);
    }
    err_control.get_errD() =Kjm2->norm()/sqrt(ode->get_system_size());
    
    n_f_eval=n_f_eval+s;
}  

// ----------------------------------------------------------------------------

RKU1::RKU1(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_,ode_)
{
    Astable = false;
    err_order = 1;
    
    damping = 0.0;
    
    reinit_integrator();
}

RKU1::~RKU1()
{

}

void RKU1::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by RKC. It follows the formula
     * \f$\rho \Delta t \leq \beta s^2 \f$ with some cautions.
     */
    
    unsigned int stages_limit = 1e6;
    unsigned int safe_add = 2;

    if(h*eigmax<1.5)
        s = 1;
    else if(h*eigmax<2.) //for s=1 we have beta=2 and not 2-4/3*damping, so maybe one stage is enough
        s=1+safe_add;
    else
    {
        s = ceil(sqrt(1+1.5*h*eigmax)-1.);
        if(damping>0.)
            while(h*eigmax >= ChebyshevMethods::ls_RKU1(s,damping))
                ++s;
        s += safe_add;
    }

    //we try to avoid more than max_s stages for roundoff errors reasons
    if(s>stages_limit)
    {
        //if more than max_s stages are needed we resize dt
        s=stages_limit; 
        if(damping==0.)
            h=0.9*(2./3.)*s*(s+2.)/eigmax;     
        else
            h = 0.9*ChebyshevMethods::ls_RKU1(s,damping)/eigmax;
        last=false;
        err_control.update_hn(h);
    }
    else if(dt_adaptivity)
        s = max(s,2); // at least two stages, needed to estimate the error
            
    if(s>s_max)
        s_max=s;    
    s_avg += s;
}

void RKU1::step(const Real t, const Real& h)
{
    Vector*& fn= integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> d;
    vector<Real> nu;
    vector<Real> mu;
    vector<Real> kappa;
    
    ChebyshevMethods::CoefficientsRKU1(mu, nu, kappa, c, d, s, damping);
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    ode->f(t,*Kjm1,*fn); 
    *Kj = mu[0]*h*(*fn);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        ode->f(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);                                   
    }   
    
    swap_ptr=ynpu;
    ynpu=Kj;
    Kj=swap_ptr;

    //compute errors
    Real mult, a1, a2, a3;    
    if(s>2)
    {
        mult = (d[s-1]-1.)/((c[s-3]-1.)*(d[s-2]-d[s-1])-(c[s-2]-1)*(d[s-3]-d[s-1]));
        a1 = (1.-c[s-2])*mult;
        a2 = (c[s-3]-1.)*mult;
        a3 = (c[s-2]-c[s-3])*mult;        
    }
    else// if(s==2)
    {        
        mult = (d[s-1]-1.)/d[s-1]/c[s-2];
        a1 = (1.-c[s-2])*mult;
        a2 = -mult;
        a3 = c[s-2]*mult;
    }
    *Kjm2 *= a1;
    *Kjm2 += a2*(*Kjm1);
    *Kjm2 += a3*(*ynpu);
     
    for(int i=0;i<ode->get_system_size();i++)
        (*Kjm2)(i) = (*Kjm2)(i)/(a_tol+r_tol*max(abs((*ynpu)(i)),abs((*yn)(i))));
    
    err_control.get_errD() =Kjm2->norm()/sqrt(ode->get_system_size());
            
    n_f_eval=n_f_eval+s;
}  


RKU2::RKU2(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_,ode_)
{
    Astable = false;
    err_order = 2;
    
    damping = 0.;
    
    reinit_integrator();
}

RKU2::~RKU2()
{

}

void RKU2::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by RKC. It follows the formula
     * \f$\rho \Delta t \leq 0.65 (s^2-1) \f$ with some cautions.
     */
    
    unsigned int stages_limit = 1e6;
    unsigned int safe_add = 2;
    
    s = ceil(sqrt(4.+2.5*h*eigmax)-1.);
    if(damping>0.)
        while(h*eigmax >= ChebyshevMethods::ls_RKU2(s,damping))
            ++s;
    s += safe_add;
                
    //we try to avoid more than max_s stages for internal stability reasons
    if(s>stages_limit) 
    {
        //if more than max_s stages are needed we resize h       
        s = stages_limit;
        if(damping==0.)
            h = 0.9*0.4*(s-1.)*(s+3.)/eigmax;
        else
            h = 0.9*ChebyshevMethods::ls_RKU2(s,damping)/eigmax;
        last=false;
        err_control.update_hn(h);
    }
    else//at least 2 stages
        s=max(s,2);
            
    if (s>s_max)
        s_max=s;
    s_avg += s;
}

void RKU2::step(const Real t, const Real& h)
{
    /**
     * Does a RKC step. Given \f$ y_n\f$ at \f$t\f$ it computes 
     * \f$y\f$ at \f$t+h\f$. Parabolic provides the right hand side function
     * and sets the Dirichlet boundary conditions.
     */
    
    Vector*& fn= integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> mu;
    vector<Real> nu;
    vector<Real> kappa;
    vector<Real> gamma;
    
    ChebyshevMethods::CoefficientsRKU2(mu, nu, kappa, gamma, c, s, damping);
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    ode->f(t,*Kjm1,*fn); 
    *Kj = mu[0]*h*(*fn);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        ode->f(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);     
        *Kj += gamma[j-1]*h*(*fn);
        *Kj += (1.-nu[j-1]-kappa[j-1])*(*yn);
    }   
    
    swap_ptr=ynpu;
    ynpu=Kj;
    Kj=swap_ptr;

    // computing local error
    ode->f(t+h,*ynpu,*Kjm1); 
    for(int i=0;i<ode->get_system_size();i++)
    {
        mu[0] = max(abs((*ynpu)(i)),abs((*yn)(i)));
        nu[0] = 0.8*((*yn)(i)-(*ynpu)(i))+0.4*h*((*fn)(i)+(*Kjm1)(i));
        (*Kjm2)(i) = nu[0]/(a_tol+mu[0]*r_tol);
    }
    err_control.get_errD() =Kjm2->norm()/sqrt(ode->get_system_size());
    
    n_f_eval=n_f_eval+s;
}  


ROCK2::ROCK2(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_,ode_)
{
    Astable = false;
    err_order = 2;
    
    reinit_integrator();
}

ROCK2::~ROCK2()
{   
}

void ROCK2::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by ROCK2. It follows the formula
     * \f$\rho \Delta t \leq 0.811 s^2 \f$ with some cautions.
     * Remark: in the algorithm s[0] isn't the mathematical \f$s\f$, actually
     * s[0]=\f$s-2\f$.
     */
    
    unsigned int safe_add = 2;
    
    s=(int)sqrt((1.5+h*eigmax)/0.811); //here s actually is s-1  
    s = s + safe_add;        
    
    if(s>199) //s>200, s can't be higher than 200 with ROCK2. We don't have the coefficients.
    {
        s=198; //s=s-2 with s=200
        h=25950.8/eigmax; //25950.8=0.8*(200*200*0.811-1.5), here s=200
        last=false;
        err_control.update_hn(h);
    }
    else
        s=max(s,2)-1; //s is s-2
    
    if(s!=sold) //if the number of stages has changed we recompute mp, which
        mdegr(s,mp);//is used to find the right coefficients in the tables
    
    if (s+2>s_max)
        s_max=s+2;    
    s_avg += s;
    
    sold = s;
}

void ROCK2::step(const Real t, const Real& h)
{
    /**
     * Does a ROCK2 step. Given \f$ y_n\f$ at \f$t\f$ it computes 
     * \f$y\f$ at \f$t+h\f$. Parabolic provides the right hand side function
     * and Dirichlet boundary conditions.
     */
    
    Real temp1,temp2,temp3;
    Real ci1,ci2,ci3;
    int mr,mz;
       
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
//             Initializations
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
    Vector*& fn= integr[0];
    Vector*& yjm1= integr[1];
    Vector*& yjm2= integr[2];
    Vector*& yj= integr[3];
    Vector* swap_ptr=0;
 
    mz=mp[0];//used to look in tables
    mr=mp[1];//as well
    
    temp1=h*recf[mr];   //temp1=dt*mu_1
    ci1=t+temp1;        //ci1 = t+dt*mu_1 = t+dt*c_1
    ci2=ci1;            //ci2 = t+dt*mu_1
    ci3=t;              //ci3 = t
    
    
    //initialization
    ode->f(t,*yn,*fn);
    
    //first step
    if(s>=2)
    {
        *yjm2=*yn;
        *yjm1=*yn+temp1*(*fn);
    }
    else
    {   //in this case we do not enter the next loop
        *yj=*yn+temp1*(*fn);    
    }

    //Stages for j=2...s-2 (remember that for rock2 s is s-2)
    for(int j=2;j<=s;j++)
    {
        temp1= h*recf[++mr];     //temp1 = dt*mu_j
        temp3= -recf[++mr];      //temp3 = -kappa_j
        temp2= 1.0-temp3;        //temp2 = -nu_j
        ci1=temp1+temp2*ci2+temp3*ci3; //ci1 = t+dt*c_j

        //computing yj=g_j
        ode->f(ci2,*yjm1,*yj);    //yj = f(g_{j-1)})
        *yj *= temp1;            
        *yj += temp2*(*yjm1);
        *yj += temp3*(*yjm2);

        //Shift the value "y" for the next stage
        if(j<s)
        {
            swap_ptr=yjm2;
            yjm2 = yjm1;   //yjm2 = g_{j-1}
            yjm1 = yj;      //yjm1 = g_j
            yj = swap_ptr;  //yj = free space previously occupied by yjm2
            
            ci3=ci2;        //ci3 = t+dt*c_{j-1}
            ci2=ci1;        //ci2 = t+dt*c_j
        }
    }
    // at this point we have yj = g_{s-2}, ci1 = t+dt*c_{s-2}
    // and yjm1, yjm2 are usable working spaces
    
    //Begin of the two-stage finishing procedure
    temp1=h*fp1[mz];     //temp1 = dt*sigma
    temp2=h*fp2[mz];     //temp2 = -dt*sigma*(1-tau/sigma^2)
    
    ode->f(ci1,*yj,*yjm2);//yjm2 = f(g_{s-2})

    ci1=ci1+temp1;        //ci1 = t+dt*c_{s-1}
    *yj += temp1*(*yjm2);  //yj = g_{s-1}

    ode->f(ci1,*yj,*yjm1); //yjm1 = f(g_{s-1})
      
    ci1 = ci1+temp1;      //ci1 = t+dt*cs = t+dt
    
    //Next we compute the two last stages and the local error
    *yjm2 *= -temp2;
    *yjm2 += temp2*(*yjm1); //yjm2 = -dt*sigma*(1-tau/sigma^2)*(f(g_{s-1})-f(g_{s-2})
    *yj += temp1*(*yjm1);   //y = g_s^*
    *yj += *yjm2;           //y = g_s

    //compute errors
    for(int i=0;i<ode->get_system_size();i++)
    {
        ci1 = max(abs((*yj)(i)),abs((*yn)(i)));
        (*yjm1)(i) = (*yjm2)(i)/(a_tol+ci1*r_tol);
    }
    (*yjm1)/= sqrt(ode->get_system_size());
    err_control.get_errD() = yjm1->norm();
    
    swap_ptr = ynpu;
    ynpu=yj;
    yj=swap_ptr;
    
    //update the number of right hand side evaluations
    n_f_eval=n_f_eval+s+2;
}  

void ROCK2::mdegr(int& mdeg, int mp[])
{
    /**
     * Given mdeg=s-2, where s is the optimal stages number, this function 
     * gives the smallest available number of stages which is greater than s.
     * (we don't have coefficients sets for every stage) 
     * The the mp variable tells where the coefficients set is located in the
     * recf array.
     */

    mp[1]=0;

    for(int i=1;i<=46;i++)
    {
        if(ms[i-1]>=mdeg) //searching for smallest compatible s
        {
          mdeg=ms[i-1];
          mp[0]=i-1;
          break;
        }
        mp[1]=mp[1]+ms[i-1]*2-1;//jump to next coefficients set
    }
}


DROCK2::DROCK2(Parameters* param_, Ode* ode_)
:ROCK2(param_,ode_)
{
}

DROCK2::~DROCK2()
{
}

void DROCK2::update_n_stages_and_h(Real& h)
{
    s=(int)sqrt((1.5+h*eigmax)/0.42)-1; //here s=s-1  
            
    if(s>199) //s>200
    {
        s=198; //s=s-2 with s=200
        h=13709/eigmax; //13709=0.8*(202*202*0.42-1.5), here s=200
        err_control.update_hn(h);
        last=false;
    }
    else //minimum 4 stages
        s=max(s,3)-1; //s[0]=s-2
        
//    cout<<"s modified"<<endl;
//    s = 11;
    
    if(s!=sold) //if the number of stages has changes we recompute mp, which
        mdegr(s,mp);//is used to find the right coefficients in the tables
    
    if (s+2>s_max)
        s_max=s+2;
    s_avg += s;
    
    sold=s;
}

void DROCK2::step(const Real t, const Real& h)
{
    /**
     * Does a damped ROCK2 step. Given \f$ y_n\f$ at \f$t\f$ it computes 
     * \f$y\f$ at \f$t+h\f$. Parabolic provides the right hand side function
     * and Dirichlet boundary conditions.
     */
    
    // NOT OPTIMAL IMPLEMENTATION, JUST TO CHECK BEFORE DOING STOCHASTICS
    
    Real temp1,temp2,temp3;
    Real ci1,ci2,ci3;
    int mr,mz;
       
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
//             Initializations
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
    Vector*& fn= integr[0];
    Vector*& yjm1= integr[1];
    Vector*& yjm2= integr[2];
    Vector*& yj= integr[3];
    Vector*& ytmp= integr[4];
    Vector* swap_ptr=0;
 
    mz=mp[0];//used to look in tables
    mr=mp[1];//as well
    const Real alpha = recalph[mz];
    const Real ha = alpha*h;
    
    
    temp1=ha*recf[mr];   //temp1=ha*mu_1
    ci1=t+temp1;        //ci1 = t+ha*mu_1 = t+ha*c_1
    ci2=ci1;            //ci2 = t+ha*mu_1
    ci3=t;              //ci3 = t
    
    
    //initialization
    ode->f(t,*yn,*fn);
    
    //first step
    if(s>=2)
    {
        *yjm2=*yn;
        *yjm1=*yn+temp1*(*fn);
    }
    else
    {   //in this case we do not enter the next loop
        *yj=*yn+temp1*(*fn);    
    }

    //Stages for j=2...s (remember that for rock2 s is s-2)
    for(int j=2;j<=s+2;j++)
    {
        if(j<=s)
        {
            temp1= ha*recf[++mr];     //temp1 = ha*mu_j
            temp3= -recf[++mr];      //temp3 = -kappa_j
            temp2= 1.0-temp3;        //temp2 = -nu_j
            ci1=temp1+temp2*ci2+temp3*ci3; //ci1 = t+ha*c_j
        }
        else// j=s+1 or j=s+2
        {
            mr = mz*4 + 2*(j-s-1);
            temp1= ha*recf2[mr];     //temp1 = ha*mu_j
            temp3= -recf2[++mr];      //temp3 = -kappa_j
            temp2= 1.0-temp3;        //temp2 = -nu_j
            ci1=temp1+temp2*ci2+temp3*ci3; //ci1 = t+ha*c_j
        }

        //computing yj=g_j
        ode->f(ci2,*yjm1,*yj);    //yj = f(g_{j-1)})
        *yj *= temp1;            
        *yj += temp2*(*yjm1);
        *yj += temp3*(*yjm2);

        //Shift the value "y" for the next stage
        if(j<s+2)
        {
            swap_ptr=yjm2;
            yjm2 = yjm1;   //yjm2 = g_{j-1}
            yjm1 = yj;      //yjm1 = g_j
            yj = swap_ptr;  //yj = free space previously occupied by yjm2
            
            ci3=ci2;        //ci3 = t+dt*c_{j-1}
            ci2=ci1;        //ci2 = t+dt*c_j
        }
    }
    // at this point we have yj = g_{s}, ci1 = t+dt*c_{s}
    // yjm1 = g_{s-1}, ci2 = t+dt*c_{s-1}
    // yjm2 = g_{s-2}, ci3 = t+dt*c_{s-2}

    //Begin of the two-stage finishing procedure
    Real sigma, tau;
    sigma=fp1[mz]; 
    temp2=fp2[mz];  //temp2 = -sigma*(1-tau/sigma^2)
    tau = sigma*(temp2+sigma);
    
    Real sigma_a, tau_a;
    sigma_a = (1.-alpha)/2.+alpha*sigma;
    tau_a = (alpha-1.)*(alpha-1.)/2.+2.*alpha*(1.-alpha)*sigma+alpha*alpha*tau;
    
    ode->f(ci3,*yjm2,*ytmp);//ytmp = f(g_{s-2})
    *yjm1 = *yjm2 + 2.*tau_a*h*(*ytmp); // yjm1 = K_{s-1}^**
    ci2 = ci3+2.*tau_a*h;
    
    ode->f(ci2,*yjm1,*yj); // yj = f(K_{s-1}^**)
    
    //Next we compute the last stage and the local error
    *yjm1 = *ytmp-*yj; // error not scaled yet
    *yj *= 0.5*h;
    *yj += *yjm2;
    *yj += (2.*sigma_a-0.5)*h*(*ytmp);
    // yj is final solution

    //compute approximate cheap local error
    for(int i=0;i<ode->get_system_size();i++)
    {
        ci1 = max(abs((*yj)(i)),abs((*yn)(i)));
        (*yjm2)(i) = (*yjm1)(i)/(a_tol+ci1*r_tol);
    }

    err_control.get_errD() = 0.5*h*(1.-sigma_a*sigma_a/tau_a)*(yjm2->norm())/sqrt(ode->get_system_size());
//    err_control->get_errD() = 0.5*h*(1.-sigma_a*sigma_a/tau_a)*norm(*yjm2)/ode->system_size();            
    
    //updating solution
    swap_ptr = ynpu;
    ynpu=yj;
    yj=swap_ptr;
    
    //update the number of right hand side evaluations
    n_f_eval=n_f_eval+s+4;
}  


// ROCK2 precomputed coefficients ----------------------------------------------

//the available stages
int ROCK2::ms[46]={1,2,3,4,5,6,7,8,9,10,
11,12,13,14,15,16,17,18,19,20,
22,24,26,28,30,33,36,39,43,47,
51,56,61,66,72,78,85,93,102,112,
123,135,148,163,180,198};

//mu and nu
Real ROCK2::recf[4476]={0.1794612899156781,0.09326607661089206,0.1268473641290642,0.02103378190528467,0.0578662712551911,0.07776456841673993,0.02157791817707098,0.09528922876588625,0.08723622960881584,0.03957407527591189,
0.05338681630900877,0.02247986572029358,0.06393794866013652,0.0861675960246663,0.07433164225960258,0.1608966052260327,0.02881884349012583,0.03903194468364506,0.02315538732440072,0.04666180548151659,
0.08821448969005685,0.05288869100008489,0.1558852247558841,0.05955157767537859,0.2309801425356076,0.02193602641069006,0.02980001095383432,0.02364036352805737,0.03566742298660092,0.0901769635289683,
0.04023410528284037,0.1575736612169138,0.0442110049482293,0.221743062150858,0.04869809870617736,0.2940499060453847,0.01726461340854506,0.02350655101084081,0.0239919520224083,0.02817607209777294,
0.09173675182820953,0.03176838133583184,0.1600219376374877,0.03469194451695863,0.2221486436776174,0.0373876350116119,0.2814712676591709,0.04053195010515956,0.349859771934868,0.01394572558988917,
0.01901965402711647,0.02425291169233543,0.02282720812165174,0.09294575674653435,0.02575013323611407,0.1622227426315588,0.0280770541578923,0.2243919049772138,0.03005377306193452,0.2802435651401261,
0.03196472247412495,0.3348331353768858,0.03423680369928859,0.398901664156519,0.01150212461258425,0.01570710434739234,0.02445103607680476,0.01887147450684578,0.09388711032641101,0.02130209293934879,
0.1640457400037581,0.02322314183593871,0.2267657098371261,0.02480419563989607,0.2818605427391582,0.02620108024481186,0.3319747646481704,0.02760343567567814,0.3822711818317003,0.02928725018247521,
0.4419889482818066,0.009650189035331108,0.01319127118048348,0.02460463437251803,0.01586242507387037,0.0946291469184127,0.01791697874611395,0.1655306545097223,0.01953741662819831,0.2288845437800309,
0.0208544873661066,0.2840687821706643,0.02197276738303543,0.332762163420867,0.02299581340430604,0.3779528907207973,0.02405334751606276,0.4244163000944219,0.02532773055505536,0.4799238174596328,
0.008212831461177801,0.01123533801931148,0.0247259253919068,0.01351988263110243,0.09522208667469109,0.01527968649293853,0.1667413647253854,0.01666740304836073,0.2306916351605186,0.01778972104750926,
0.2862271893660237,0.01872598244401872,0.3345934460509904,0.01954407472805166,0.3778520013827246,0.02031544030249,0.4188775043096391,0.02113066217102199,0.4619185110779211,0.02211266742157054,
0.5134449754977242,0.007074719613815112,0.009684509591053977,0.02482326781574381,0.01166040318768419,0.09570219065303219,0.0131845989014643,0.1677351729038757,0.01438717215812115,0.2322139970245483,
0.01535795943897448,0.2881631572191675,0.01616136143489047,0.3366119552764601,0.01684787625840699,0.3791920530521655,0.01746359793359561,0.4179163657739238,0.01805925334195525,0.4554086296687285,
0.01869919690078058,0.495380320313383,0.01946782991114719,0.5431953913817005,0.006158116280205625,0.008434133870779752,0.02490250308471492,0.01015973796908432,0.09609568511780181,0.01149255448537381,
0.1685578212579456,0.01254496388900155,0.2334953976202505,0.0133941202672405,0.2898496927042313,0.01409430683259744,0.3385311597112782,0.01468607650095395,0.380955385363272,0.01520306181269769,
0.4187183922578105,0.01567773486210682,0.4536627812581018,0.01614693417033556,0.4881324072137513,0.01665708022899876,0.5253402953586265,0.0172671276878042,0.569731047549543,0.005408899288459016,
0.007411165652735432,0.02496784203339785,0.008931010375660927,0.09642194989632852,0.01010622203571924,0.1692450860120974,0.01103497280658699,0.2345784420179126,0.01178444954979593,0.2913056300783272,
0.01240145878830343,0.3402656182942943,0.01292005530031986,0.382756777938251,0.0133668781100104,0.4201602339064523,0.01376523513966666,0.4539235794662016,0.01413868046925951,0.4856954694207845,
0.01451444200099929,0.5175484724211766,0.01492644572980292,0.5522503690808819,0.01541651122328756,0.5934874290489524,0.004788707289765767,0.006563722805402352,0.0250223058569931,0.007912416503013329,
0.09669512238638367,0.008956313102014822,0.1698239427218661,0.009781969024648538,0.2354985079518828,0.01044851744036037,0.2925601105800896,0.01099693219507985,0.3418012690877302,0.0114565920869984,
0.3844521174650939,0.01184970544276362,0.4217763245489749,0.01219446847463953,0.4550117593290032,0.01250758855429112,0.4854408050495134,0.01280655793511545,0.514526804170333,0.01311180392585363,
0.54409412649029,0.01344843304771755,0.5765229845165943,0.01384655622092475,0.6148795064211108,0.004269409473142479,0.005853691221954542,0.02506818803768801,0.007058488452666192,0.09692609816181726,
0.007991811416685026,0.170315738084356,0.008730579305889177,0.2362853743226745,0.009327287399213484,0.2936438651173445,0.009818197983040462,0.3431514772416978,0.01022909697576184,0.3859958359229893,
0.01057909230023007,0.4233725763238417,0.01088320687501712,0.4564025740290915,0.0111543058671517,0.4861671131891395,0.01140470198717317,0.5137949407970326,0.01164762828827466,0.5405812657784165,
0.01189859698528282,0.568127910648508,0.0121763907812416,0.5984802463148428,0.0125029779529791,0.6341959825388761,0.003830262680422776,0.00525292423946998,0.02510718258012232,0.006335607341496969,
0.09712300411103621,0.007174964200736387,0.1707366515477548,0.007839821469171155,0.2369623770049865,0.008377122084280812,0.2945833847967025,0.008819233466500509,0.3443363340392364,0.009189060633588098,
0.3873807210303031,0.00950337880056659,0.4248706039873364,0.00977504705184527,0.4578585641947919,0.01001457306569686,0.4873111192988822,0.01023132837651919,0.5141675256404127,0.01043459649534057,
0.5394214279393358,0.01063454146082206,0.5642188354052674,0.01084307018096031,0.5899649430016252,0.01107438039614643,0.6184184090843611,0.01134470195438803,0.6517206797990772,0.003454097785831183,
0.004738087690909603,0.02514098616718916,0.005715864568881495,0.09729414530950285,0.006474399737140852,0.1711037123670141,0.007075625560198937,0.2375553213635903,0.007561765789488958,0.2954111541648083,
0.007961900378019847,0.3453896252983967,0.008296546652001753,0.388630302094808,0.008580623356383654,0.4262600339398994,0.008825386359541219,0.459288916256125,0.009039751450513937,0.4886139630455138,
0.009231266796183848,0.5150610315193943,0.009406897826583457,0.5394434170928742,0.009573720059444631,0.5626311157802961,0.009739556257960885,0.5856281920506604,0.009913518363761172,0.6096518917386764,
0.01010629643624125,0.6361953557375821,0.01032986385558493,0.6670333276830246,0.003131660281084644,0.004296601495430642,0.02516987418923732,0.005184210888756367,0.09744072529384443,0.005873198129887241,
0.1714189769006723,0.006419613566950146,0.2380664100570842,0.006861660893085171,0.2961280333876439,0.007225632384805198,0.346308042021841,0.007530044942257158,0.3897316580976775,0.007788304877568383,
0.4275077537500284,0.008010430934156831,0.4606200718057599,0.008204207591967733,0.4899241115804965,0.008376002129523446,0.5161793704075401,0.008531389917004646,0.5400945718660145,0.008675676859844068,
0.5623795630774431,0.008814368379677447,0.5838024839647016,0.0089535958165003,0.6052506199023558,0.009100459154569437,0.6277890216584121,0.009263166205736671,0.6527013571399886,0.00945074055537075,
0.6814810379908328,0.002843902743022214,0.003902467368190393,0.02519598783827689,0.004709429895823432,0.09757348323445655,0.005336152661315786,0.1717052034351406,0.005833455676281756,0.2385318286257855,
0.006235971415149352,0.2967834040784211,0.006567522559559145,0.3471521575495951,0.006844866528765064,0.3907519751170857,0.007080101687099548,0.4286785748173565,0.007282213357366904,0.461897608939992,
0.007458094741776693,0.4912372216070018,0.007613251892383208,0.5174141871126332,0.007752321382046524,0.5410706897885087,0.007879480568321315,0.5628149979263304,0.007998799282231232,0.5832647350414633,
0.008114557939046675,0.6030926481771959,0.008231531848213921,0.623073488760433,0.008355206753895566,0.6441267697250795,0.008491839744291131,0.6673428282976036,0.008648214015485161,0.6939680907913989,
0.002388386682754135,0.003278278227000989,0.02523707236937346,0.003957198838737114,0.09778285117940604,0.004484932501394734,0.1721579526332683,0.004904057483051761,0.2392707453970052,0.005243583901636578,
0.2978287368553751,0.00552345600742706,0.348506768678122,0.005757685183723855,0.3924033610212076,0.005956358204781419,0.4305980352749587,0.006126919713095007,0.4640361265644109,0.006275007901630985,
0.4935170207127558,0.006405016526854853,0.5197143572851119,0.00652048877918682,0.5432046533152863,0.006624408123666273,0.5644970257796814,0.006719427053538581,0.5840623664365776,0.006808059512718916,
0.6023621974962409,0.006892851836922015,0.6198777517486085,0.006976537015760206,0.6371392310474401,0.007062164945929763,0.6547536993606798,0.007153184521396811,0.6734273362784953,0.007253430606435754,
0.6939734212185117,0.007366943161281851,0.7172914790760468,0.002034217998243191,0.002792734020582961,0.02526908982706098,0.003371787305020115,0.097946430626842,0.00382220106251195,0.1725128108619325,
0.004180174912236582,0.2398521447789659,0.004470370456915701,0.2986551495007976,0.00470973753959748,0.3495841277608567,0.004910176531706057,0.393727215546069,0.005080243210917943,0.4321540636611387,
0.005226234784652425,0.4657990145921906,0.005352893648837146,0.4954473934553277,0.005463875356202293,0.521753178861748,0.005562069670587001,0.5452642338468077,0.005649829175907964,0.5664474938420431,
0.005729139522573558,0.5857122023051498,0.005801753106324432,0.6034312542567039,0.005869300229750696,0.6199612865321288,0.005933386406539511,0.6356621535317655,0.005995679940966749,0.650916092177813,
0.006057989063542723,0.6661462081376089,0.006122321685397079,0.681832779782845,0.00619091247673213,0.6985240974244095,0.006266191569573649,0.7168360261395987,0.006350658968555684,0.7374314365934664,
0.001753390300533029,0.002407592405123618,0.025294539852267,0.002907262330175829,0.09807671703676038,0.00329614305335325,0.172796146245125,0.003605393786629613,0.240317752902561,0.003856241133152162,
0.2993193629853389,0.004063271707574896,0.3504538499243808,0.004236723270781009,0.3948018542411614,0.004383951358482999,0.433426391661283,0.004510362519918408,0.4672550951343391,0.00462001758635921,
0.4970655128750277,0.004716030788908096,0.5235015776135139,0.004800840911719355,0.5470972903869747,0.004876401000932111,0.5682998788764699,0.004944315587632199,0.5874903820382004,0.005005943872710219,
0.6050015924522462,0.005062480886612362,0.6211339094140743,0.005115024538196171,0.6361697673414449,0.005164633649628788,0.6503872077614197,0.005212379842880342,0.664072947534986,0.005259393975018879,
0.6775349550821846,0.005306905271212799,0.6911140242391297,0.005356268008877997,0.7051930501274516,0.005408966392479879,0.7202015894622125,0.005466583459521371,0.73661183744624,0.005530715851949468,
0.7549206082556241,0.001526985597447617,0.00209699647376184,0.02531509673360445,0.002532539440756325,0.09818212210216928,0.002871665019548934,0.1730258242612469,0.00314147957617889,0.2406960805778798,
0.003360448272795954,0.2998605891705305,0.003541258850863859,0.351164920618253,0.003692814597959592,0.3956840637564039,0.00382150991829634,0.4344762665775814,0.003932041593754749,0.4684646795111045,
0.004027932600409154,0.4984220865031195,0.00411187783441815,0.5249867454711314,0.004185977901331621,0.5486853051737176,0.004251901290680174,0.5699549871771356,0.004310999960360177,0.5891629130044663,
0.004364394207792517,0.6066224358437886,0.004413037135908646,0.6226069634301489,0.004457765549253931,0.6373618910220726,0.004499341872619412,0.6511152177198366,0.004538490146036315,0.6640873078294923,
0.004575927973865173,0.6765001110571239,0.004612395244609121,0.6885859633887566,0.004648679290544145,0.7005958269983906,0.004685634750981115,0.7128064532627977,0.004724194620734037,0.725525427370265,
0.004765366815800576,0.7390923556807402,0.004810208361466176,0.7538736313498631,0.00485976779803573,0.7702474410767316,0.001341778429126211,0.001842856050972565,0.02533193854590164,0.002225855496752996,
0.09826859006407132,0.002524180832845478,0.173214539260328,0.002761629971248482,0.2410075282122661,0.002954414092304398,0.3003071440620802,0.003113671483601745,0.3517531700231469,0.003247217194092525,
0.3964161974860505,0.00336066282214273,0.4353508843394245,0.003458128531982532,0.4694771791546352,0.003542702530515273,0.4995646943102908,0.003616743969933074,0.5262481851213445,0.003682087306007793,
0.5500500545141821,0.003740183451080423,0.5714020312900167,0.003792199619292558,0.5906638632386166,0.003839091722399387,0.6081388453953844,0.003881658284331042,0.6240866310409383,0.003920581807177477,
0.6387339085746301,0.003956461592009382,0.6522834922619609,0.003989840752967416,0.6649222876852415,0.004021229291927125,0.6768284931080795,0.004051124448454234,0.6881782931898645,0.0040800289816473,
0.6991521843388895,0.004108467479283411,0.7099409264259735,0.004137000143073717,0.720750923666178,0.004166232701476242,0.7318085765567045,0.004196820124371787,0.7433628001322098,0.004229460707743591,
0.7556844737332599,0.004264876064174216,0.7690611204503529,0.004303772042593957,0.7837847386079905,0.001121314683543782,0.001540263996358294,0.02535201582246297,0.001860616782749459,0.09837180081070002,
0.002110257847461394,0.1734401492675452,0.00230905472549513,0.2413805649139678,0.002470541128098822,0.3008431837086572,0.002604015187479208,0.3524611103987945,0.002716000109093522,0.3972999402435936,
0.00281117899964196,0.4364103598053284,0.002892989192294541,0.4707089285552096,0.002964005631609554,0.5009621096861159,0.003026193412069002,0.5278014394899639,0.003081077890619739,0.5517457095527168,
0.003129861843706196,0.5732222749918362,0.003173507910562157,0.5925852995229407,0.003212797846468413,0.6101307283827272,0.00324837602242407,0.6261084073494506,0.003280782072682338,0.6407319017471913,
0.003310475987498136,0.6541865362303725,0.00333785791196733,0.6666360953708608,0.003363284226199348,0.6782285398626227,0.003387081014309363,0.6891010163282705,0.003409555695216232,0.6993843710724037,
0.003431007330102873,0.7092073156356478,0.003451735898380667,0.7187003283870178,0.003472050613465582,0.7279993040305766,0.003492277104239913,0.7372488732478574,0.003512762995683089,0.7466051990304045,
0.003533881070388274,0.7562379075535515,0.003556028788751684,0.7663306282943305,0.003579622533335547,0.7770794116956159,0.003605084625764653,0.7886880974410957,0.003632821128406038,0.8013595930224617,
0.0009510596127542824,0.001306530384822369,0.02536754145655561,0.001578428318009605,0.09845171104183326,0.001790384464545768,0.1736150904003335,0.001959236797520781,0.2416703438788951,0.002096455322227467,
0.3012604613663233,0.002209920074417597,0.3530135437139923,0.002305159021830049,0.3979914926982783,0.002386140355447109,0.4372421218440167,0.002455776174573485,0.4716796241323911,0.002516247355295036,
0.5020683797660233,0.002569217475366734,0.5290379452961959,0.002615976824974766,0.5531050686294569,0.00265754147041669,0.5746948225947125,0.002694722818617351,0.594158651817051,0.002728177434091533,
0.6117891098846177,0.002758443393525188,0.6278316916720729,0.002785967312803683,0.6424943013679313,0.002811124820961013,0.6559548636293389,0.002834236378150614,0.6683675049794192,0.002855579757979202,
0.6798676490616872,0.002875400127857228,0.6905762965442401,0.002893918395849096,0.7006037006910071,0.002911338305155523,0.7100526010841707,0.002927852618951001,0.7190211375279999,0.002943648628770056,
0.7276055302511653,0.002958913124102285,0.7359025774945679,0.00297383686744952,0.7440119836573593,0.002988618517881438,0.7520384865047949,0.003003467829003067,0.7600936967329708,0.003018607809042993,
0.7682974943112553,0.003034275371809675,0.7767787422315526,0.003050719837809893,0.785674982357174,0.003068198491381527,0.7951306800776951,0.003086968312155531,0.8052935061547398,0.00310727305576029,
0.8163081231039382,0.0008168405299925574,0.001122235164426433,0.02537979319269659,0.001355887567479481,0.09851483086228864,0.001538080541563495,0.1737534362414166,0.001683267147792567,0.2418998259562802,
0.001801292850812278,0.3015914511597799,0.001898921476480742,0.3534525616677698,0.001980897763334222,0.3985422441814719,0.00205062745022216,0.4379061538920372,0.002110609627253456,0.4724567457971343,
0.002162715579344156,0.5029569229728563,0.00220837237888526,0.5300349067999921,0.002248686443403842,0.5542061521610998,0.002284528486070424,0.5758943864439198,0.002316593115008434,0.5954495572468081,
0.002345441446771229,0.6131624579583983,0.002371532122833287,0.6292764292152166,0.002395244271202193,0.6439966702794445,0.002416894786788224,0.6574976613071741,0.002436751550587987,0.6699291171394233,
0.002455043712749272,0.6814208098492534,0.002471969833678045,0.6920865249452547,0.002487704452255031,0.7020273575786585,0.002502403494283208,0.7113345090150194,0.002516208823898241,0.7200917075437527,
0.002529252160257252,0.7283773493344921,0.002541658520726186,0.736266431138189,0.002553549302394987,0.7438323260356781,0.002565045070295341,0.7511484336838515,0.002576268078538789,0.7582897158118264,
0.002587344505847411,0.7653341042347025,0.00259840633642072,0.772363740630846,0.002609592758562768,0.7794659733020046,0.002621050886496552,0.7867339953288847,0.002632935537833889,0.7942669615898086,
0.002645407727623461,0.8021693721981406,0.002658631484279102,0.8105494642340256,0.002672768576741324,0.8195163250963701,0.00268797079956683,0.8291754492311013,0.0006779619032019482,0.0009315109162316307,
0.02539248163813982,0.001125546650705646,0.0985802566139559,0.001276892135791717,0.1738969874989156,0.001397535825902466,0.2421382412064102,0.001495644371531161,0.3019358267641449,0.001576828175558173,
0.3539100964089328,0.001645022736116236,0.3991173116956755,0.001703052777188943,0.4386009869035954,0.001752990885928625,0.4732718824089663,0.001796388809395077,0.5038914977622926,0.001834429731186084,
0.5310868327976426,0.001868030737964622,0.5553722186397624,0.001897913248880432,0.5771702844407484,0.001924652401977583,0.596829836122438,0.001948712334483372,0.6146404100545767,0.001970471824369135,
0.6308438946964408,0.001990243228225435,0.6456437499005518,0.002008286680798674,0.659212320773459,0.002024820896049421,0.6716966627198852,0.002040031498715453,0.6832232108862492,0.00205407754078053,
0.6939015548982262,0.002067096670752565,0.703827521382314,0.00207920929501315,0.7130857211033993,0.002090521980414542,0.7217516823653617,0.002101130283221459,0.7298936653027928,0.002111121143106236,
0.7375742308368269,0.002120574946574002,0.744851621733151,0.002129567338002603,0.7517810000528958,0.002138170835613056,0.7584155742314791,0.002146456292026427,0.7648076391248479,0.002154494222946041,
0.7710095427995872,0.002162356011555635,0.7770745838488886,0.002170114979286494,0.7830578318730199,0.002177847294739602,0.7890168508176265,0.002185632671118797,0.7950122895929717,0.002193554778405604,
0.801108286542413,0.002201701270375675,0.8073726141153927,0.002210163300386258,0.8138764685711868,0.002219034377453016,0.8206937889887616,0.002228408401671311,0.8278999743285053,0.002238376724427438,
0.8355698630075664,0.002249024115421597,0.8437748549346639,0.0005717171765067005,0.0007855819755678591,0.02540219552841172,0.0009492803029534314,0.09863038311884588,0.001076991591417719,0.1740070746137774,
0.001178820912278011,0.2423212835361945,0.001261652009696089,0.3022005646620261,0.001330213911048136,0.3542623490383656,0.001387823796651901,0.3995607968545022,0.001436862514566579,0.4391378457027524,
0.001479076886266864,0.4739030248719332,0.0015157746727935,0.5046168340722528,0.001547952929687645,0.5319054272436042,0.001576384371750672,0.5562823875920914,0.001601676728454795,0.5781696502895681,
0.001624314356699193,0.5979153431712045,0.001644687957566232,0.6158083043507766,0.001663116161966728,0.6320896668598296,0.001679861458243818,0.6469620376959935,0.001695142117261203,0.6605967661879515,
0.001709141243106296,0.6731397164733577,0.001722013731055875,0.6847158755342792,0.001733891682879084,0.695433055919974,0.001744888672268816,0.7053848938737558,0.001755103144742999,0.7146532979347533,
0.001764621160512043,0.7233104679728801,0.001773518635044922,0.7314205777797613,0.001781863193450182,0.7390411938444027,0.001789715726676655,0.7462244872361469,0.001797131716770826,0.7530182833863068,
0.001804162382812431,0.7594669850655551,0.001810855687135456,0.7656123962667175,0.001817257231900906,0.7714944684438817,0.001823411068174413,0.7771519851570839,0.001829360432715894,0.7826231962275617,
0.001835148421173907,0.7879464076688663,0.001840818599844342,0.7931605286026111,0.001846415551221751,0.7983055707993119,0.001851985340937259,0.8034230901440029,0.001857575885138658,0.8085565520117598,
0.001863237187886467,0.8137515941551513,0.001869021407926604,0.8190561513413812,0.001874982703844713,0.8245203960095123,0.001881176797227895,0.8301964394574721,0.001887660186864975,0.8361377299185951,
0.001894488945824785,0.8423980794999191,0.001901717040890795,0.8490302543125972,0.001909394134299323,0.8560840749607826,0.0004886273679748824,0.000671443795464978,0.02540979640585205,0.0008113983041306847,
0.09866962915608751,0.0009206047809355559,0.1740933291741878,0.001007696615844565,0.2424648238198073,0.001078555077304152,0.302408378783803,0.001137220495111209,0.3545391775646764,0.001186526928471642,
0.3999097732746233,0.001228508433020824,0.4395609067610286,0.001264657346363223,0.4744011819765317,0.00129609086494968,0.5051903569581037,0.001323660708324793,0.5325539695623672,0.00134802690605591,
0.557005071975368,0.001369708506508536,0.5789651232823037,0.001389119124621932,0.5987818060211025,0.001406592323781286,0.6167435227777786,0.001422400048200691,0.6330909609625438,0.001436766218460017,
0.6480262517065541,0.001449876904280157,0.6617202166442876,0.001461888037961801,0.6743181164063009,0.001472931335824075,0.6859442313853145,0.001483118897071317,0.6967055330639804,0.001492546815085222,
0.7066946457906411,0.001501298043441827,0.7159922532405723,0.001509444694127518,0.7246690686840157,0.001517049899494373,0.7327874613555573,0.00152416933654122,0.7404028107621932,0.001530852488184009,
0.7475646451451914,0.001537143698611175,0.7543176083360913,0.001543083066767105,0.7607022900146448,0.001548707212186899,0.7667559471993193,0.001554049939909855,0.7725131391565733,0.00155914282537307,
0.7780062934020588,0.001564015735548671,0.7832662167717693,0.001568697298775154,0.7883225624132202,0.001573215332462839,0.7932042607791593,0.001577597234898074,0.7979399201221365,0.001581870344542036,
0.8025581994310867,0.001586062267357487,0.8070881540792818,0.001590201169665478,0.8115595515386727,0.00159431603072655,0.8160031512482756,0.001598436845587962,0.820450939020954,0.001602594764727253,
0.8249363021965112,0.001606822152724702,0.829494127135844,0.001611152543807615,0.8341607957498925,0.001615620467991531,0.8389740528796195,0.001620261118280227,0.8439727120207893,0.001625109827821845,
0.849196163928842,0.001630201327174513,0.8546836521716283,0.001635568757273152,0.8604732831634605,0.001641242424784281,0.8666007472872529,0.0004079721290191449,0.0005606390646000299,0.02541717774246925,
0.0006775302011759069,0.09870776119617666,0.0007687559860220674,0.1741771883435951,0.0008415222726271052,0.2426044825025154,0.0009007379123407699,0.3026107489587744,0.00094977519824073,0.3548090209828054,
0.0009909995793709382,0.4002503224647114,0.001026108754684635,0.4399742629913482,0.001056348256439421,0.474888578965461,0.001082650488031221,0.5057523403210779,0.001105726258156915,0.5331905211809378,
0.001126126368448291,0.5577156967903346,0.001144283937532252,0.5797489091761477,0.001160544069671673,0.5996374632935472,0.001175185037186027,0.6176694070801912,0.001188433661215005,0.6340850818085325,
0.001200476654074482,0.6490862674574625,0.001211469103358439,0.6628434159685319,0.001221540901759423,0.6755013855037041,0.001230801679415364,0.6871840056733981,0.001239344630375328,0.6979977314932807,
0.001247249512545108,0.7080345854591737,0.001254585023077141,0.717374541493969,0.001261410697040181,0.7260874694062653,0.001267778438847666,0.7342347316708526,0.001273733768408534,0.7418705038821782,
0.001279316843997812,0.7490428746151535,0.001284563309197927,0.7557947684726961,0.001289505000406732,0.7621647269098429,0.001294170543282941,0.7681875743244604,0.001298585860357731,0.7738949913863433,
0.001302774607351322,0.7793160132569095,0.001306758552112867,0.7844774669407445,0.001310557907271291,0.7894043592864977,0.001314191625436819,0.7941202249469652,0.001317677663971441,0.7986474417822638,
0.001321033224831176,0.8030075196380616,0.001324274973678666,0.8072213670641986,0.001327419241294216,0.8113095392825265,0.001330482209212664,0.8152924695008139,0.001333480080427699,0.819190684443498,
0.001336429234887632,0.8230250036771799,0.001339346368318267,0.8268167209024656,0.001342248611619848,0.8305877638256554,0.001345153626679458,0.8343608274878259,0.001348079672919514,0.8381594740079386,
0.001351045637295423,0.8420081896122813,0.00135407101882596,0.8459323876377299,0.001357175857202399,0.8499583440311814,0.001360380593752449,0.8541130499173916,0.001363705852281417,0.8584239643605064,
0.001367172127403134,0.8629186498951948,0.001370799369307807,0.8676242742545544,0.001374606456942832,0.8725669645681554,0.001378610556753209,0.8777710057814074,0.0003457584334322357,0.0004751621453940012,
0.02542287327621289,0.0005742529802903852,0.09873719754448176,0.0006515970827455522,0.1742419598034706,0.0007132998513406398,0.2427124231693788,0.0007635205129753468,0.3027672771691929,0.0008051163841037499,
0.3550179179580828,0.0008400916805185935,0.4005142108161735,0.000869884827666364,0.4402949133244795,0.0008955511242402197,0.4752671131865228,0.0009178805496049455,0.5061893703422927,0.0009374753213346,
0.5336862457066851,0.0009548020787478288,0.5582699698916923,0.0009702277440292999,0.5803612891309488,0.0009840446600062943,0.6003072470093884,0.0009964885371329649,0.618395653440728,0.001007751484296832,
0.6348666260010463,0.001017991617422186,0.6499217274752063,0.001027340245776507,0.6637311918729225,0.001035907317136567,0.6764396515818053,0.00104378559353748,0.688170695275452,0.001051053889330778,
0.6990305140455065,0.001057779608173888,0.709110834894387,0.001064020749985821,0.7184912951109098,0.001069827513031319,0.7272413759474827,0.001075243583787733,0.7354219871919199,0.001080307183927825,
0.7430867737680106,0.001085051926828077,0.7502831998805611,0.001089507523598374,0.7570534542643864,0.001093700369430917,0.7634352109083323,0.001097654034185638,0.7694622725315966,0.001101389675935647,
0.77516511858508,0.001104926392241615,0.7805713752558724,0.001108281520887288,0.7857062215827314,0.00111147089945648,0.7905927431287384,0.001114509091294428,0.7952522425401743,0.001117409583946746,
0.7997045146224264,0.001120184965012693,0.8039680921873577,0.001122847079414237,0.8080604677961332,0.001125407171312882,0.811998295576812,0.001127876013259072,0.8157975764884198,0.001130264024599364,
0.8194738296925913,0.001132581380665879,0.8230422520462021,0.001134838113807109,0.8265178671142719,0.001137044206868379,0.8299156644956935,0.00113920967927677,0.8332507296316362,0.001141344665414018,
0.8365383636069283,0.001143459484459602,0.8397941917404019,0.001145564700346352,0.8430342589770379,0.001147671169888255,0.846275109234415,0.001149790076517388,0.8495338449183082,0.001151932946415512,
0.8528281618190179,0.001154111643169411,0.8561763535592529,0.001156338336456831,0.8595972787361906,0.001158625439740416,0.8631102829627566,0.001160985511591647,0.8667350672783336,0.001163431115191409,
0.870491494017944,0.001165974630888303,0.874399321392876,0.001168628017588916,0.8784778589727296,0.001171402520361833,0.8827455382208639,0.001174308324102304,0.8872193954732928,0.0002967633115387107,
0.0004078420646640203,0.025427359800461,0.0004929083178239419,0.09876039342539503,0.00055931267512858,0.1742930216852869,0.0006122943067439731,0.2427975603509948,0.0006554223612098552,0.3028908101007459,
0.0006911487848386307,0.3551828911240184,0.000721193448713687,0.4007227692617545,0.0007467906965300383,0.4405485438561624,0.0007688460941907184,0.47556680431992,0.0007880375248105768,0.5065357215285611,
0.0008048817439819273,0.5340795437694515,0.0008197791606092136,0.5587102435529917,0.0008330446107364092,0.5808483484024387,0.0008449289291109884,0.6008407117558368,0.0008556343498905159,0.6189749738875788,
0.0008653256884173914,0.6354910969728228,0.0008741385860890041,0.6505904974988862,0.0008821856763681789,0.6644432678330568,0.0008895612564439208,0.6771938992957999,0.000896344869325779,0.6889658361211737,
0.0009026040810144197,0.6998651175897177,0.0009083966557751288,0.7099833073204581,0.0009137722762541778,0.7193998631153847,0.0009187739158093501,0.7281840656590849,0.0009234389425261854,0.7363965975579181,
0.0009278000143753827,0.7440908437477775,0.0009318858104408874,0.7513139686822109,0.0009357216324926972,0.7581078137561323,0.0009393299032829178,0.7645096492310164,0.0009427305820375386,0.7705528078328676,
0.0009459415131589371,0.7762672216898303,0.0009489787207614479,0.781679879983595,0.0009518566590593747,0.7868152213232364,0.0009545884266143449,0.791695472197249,0.000957185950881378,0.7963409407568194,
0.0009596601482630737,0.8007702735070883,0.0009620210639090244,0.8050006811385999,0.0009642779947228965,0.8090481386453148,0.0009664395984170944,0.8129275639917815,0.0009685139909496511,0.8166529788658923,
0.0009705088342629528,0.8202376544499068,0.0009724314158979084,0.8236942446328719,0.0009742887217632853,0.8270349086491756,0.0009760875030844708,0.8302714247416513,0.0009778343383275471,0.8334152960971942,
0.000979535690683748,0.8364778499742661,0.0009811979614979323,0.8394703306225847,0.0009828275398255845,0.842403986274681,0.0009844308480999042,0.8452901501569875,0.0009860143836786285,0.8481403151159906,
0.0009875847558152445,0.8509662010755018,0.0009891487173585096,0.8537798141289356,0.0009907131902268241,0.8565934956229708,0.0009922852834316175,0.8594199591072303,0.0009938723021413953,0.8622723125149068,
0.0009954817459946895,0.8651640624148244,0.0009971212946004986,0.8681090966584658,0.0009987987779301983,0.8711216412695179,0.001000522129134342,0.8742161870354788,0.00100229931724864,0.8774073810231204,
0.001004138257331534,0.8807098782305657,0.001006046695854091,0.8841381489021722,0.001008032069698879,0.8877062377737649,0.001010101337974326,0.8914274727949142,0.001012260787062656,0.8953141227940432,
0.0002505764977838404,0.0003443770249124342,0.02543158995016821,0.0004162174537078318,0.09878227034153308,0.0004723030362070807,0.1743411977575627,0.0005170566942884504,0.2428779209209343,0.0005534915375556028,
0.3030074710503266,0.0005836775268839069,0.355338776377064,0.000609066594502194,0.4009199654813231,0.0006307007808156645,0.4407885275977928,0.0006493445644344517,0.4758505926393938,0.0006655702083843857,
0.5068639748553351,0.0006798139474548506,0.5344526379872747,0.0006924137931820147,0.5591283215877272,0.0007036355137975163,0.5813113577588189,0.0007136908453729489,0.6013484321270596,0.0007227504933228287,
0.619527037525631,0.0007309535721159746,0.6360870035914407,0.0007384145655048229,0.6512296248800408,0.0007452285316469441,0.6651248789020748,0.0007514750465678791,0.6779171461479578,0.0007572212276838003,
0.6897297612797826,0.0007625240776755495,0.7006686526246469,0.0007674323201021287,0.7108252688389807,0.0007719878506216157,0.7202789460425201,0.0007762268944500497,0.7290988336445309,0.0007801809371337034,
0.7373454702759427,0.0007838774788105678,0.7450720807914899,0.0007873406498722635,0.7523256496919639,0.0007905917169405165,0.7591478143617146,0.0007936495014053848,0.7655756123276436,0.0007965307277851307,
0.7716421096510724,0.0007992503154034614,0.7773769320588682,0.0008018216240149106,0.7828067161269425,0.0008042566618112444,0.7879554944635353,0.000806566262542953,0.7928450261873953,0.0008087602371674192,
0.7974950818949961,0.0008108475043988473,0.8019236906381128,0.0008128362037173276,0.806147355094338,0.0008147337937451987,0.8101812400362133,0.0008165471383802039,0.8140393383338362,0.0008182825826579203,
0.8177346180177749,0.0008199460199784433,0.8212791533500792,0.0008215429520572276,0.8246842423744934,0.0008230785427339746,0.8279605130215852,0.0008245576665860399,0.8311180195135344,0.0008259849531356726,
0.8341663305331678,0.0008273648273067869,0.8371146103813345,0.000828701546671456,0.8399716941366764,0.0008299992359243829,0.8427461576444749,0.0008312619189314546,0.8454463829899002,0.0008324935486128937,
0.8480806199498299,0.000833698034839669,0.8506570437613244,0.0008348792704412902,0.8531838093892227,0.0008360411553417403,0.8556691023160556,0.0008371876187562825,0.8581211857109087,0.0008383226392937222,
0.8605484436569011,0.0008394502627153142,0.8629594199271572,0.0008405746170023279,0.8653628515950123,0.000841699924279396,0.8677676965453935,0.0008428305090311761,0.8701831537222023,0.0008439708019376514,
0.8726186747045624,0.0008451253385421383,0.8750839649592925,0.0008462987518610647,0.8775889728776909,0.0008474957579532225,0.8801438644857288,0.000848721133398272,0.8827589815370717,0.0008499796836021254,
0.8854447805825755,0.0008512762008654021,0.8882117505884698,0.0008526154112376473,0.8910703067843766,0.000854001909353256,0.8940306587022231,0.0008554400807241641,0.8971026508610441,0.0008569340113668948,
0.9002955753026496,0.0008584873851806726,0.9036179562250074,0.0002143902298924797,0.0002946512237154931,0.02543490463324859,0.0003561260688438285,0.09879941719514369,0.0004041230092741482,0.1743789695511629,
0.0004424256774627034,0.2429409502831179,0.0004736116853150793,0.3030990118950347,0.0004994518169910726,0.3554611561972349,0.0005211881524129952,0.401074863083999,0.0005397121201514141,0.4409771509224894,
0.0005556777014490634,0.4760737966263269,0.0005695744414699151,0.5071223419498457,0.0005817755161210044,0.534746534570652,0.0005925700732976313,0.5594579384248554,0.0006021854581155787,0.5816767395342153,
0.0006108027921070821,0.601749499329333,0.0006185680956696988,0.6199636028202163,0.0006256003634582641,0.6365587841093024,0.0006319975186065948,0.6517362513852718,0.0006378408654613904,0.6656659024902249,
0.0006431984629598351,0.6784920429133774,0.0006481277109801403,0.6903379352397969,0.0006526773552254315,0.7013094370793644,0.0006568890572555564,0.7114979262641195,0.0006607986356285067,0.720982666549765,
0.0006644370556787454,0.7298327319935689,0.0006678312253083076,0.7381085813804109,0.0006710046397087699,0.7458633536247841,0.0006739779074395309,0.7531439394662658,0.0006767691825904827,0.7599918728233108,
0.0006793945220531763,0.7664440759826728,0.0006818681826575195,0.7725334857075471,0.000684202869710358,0.7782895818427216,0.0006864099460208997,0.7837388357016259,0.0006884996086172708,0.7889050921541919,
0.0006904810389048264,0.7938098966818099,0.0006923625308852659,0.7984727765644702,0.0006941516011688034,0.8029114836923505,0.0006958550838122305,0.807142205155724,0.000697479212460725,0.8111797466910577,
0.0006990296918282764,0.8150376931919887,0.0007005117601960531,0.8187285497884573,0.0007019302443210675,0.8222638664221132,0.0007032896079146326,0.8256543483749293,0.0007045939946601219,0.8289099548200989,
0.0007058472665836456,0.8320399871433695,0.0007070530384625938,0.8350531685159762,0.0007082147088500877,0.8379577159767848,0.0007093354882038982,0.8407614060927453,0.0007104184245328271,0.8434716351063775,
0.0007114664269090819,0.846095474341138,0.00071248228713946,0.8486397215154378,0.0007134686998392874,0.8511109485098599,0.000714428281109386,0.8535155460363737,0.0007153635859765056,0.8558597655701624,
0.0007162771247204472,0.8581497588214944,0.0007171713781754927,0.860391614944607,0.0007180488120588147,0.8625913956007935,0.0007189118903434776,0.8647551679119477,0.000719763087657754,0.8668890352571068,
0.0007206049006551845,0.8689991657766132,0.0007214398582606628,0.8710918183552704,0.0007222705306565306,0.8731733657564531,0.0007230995368291464,0.8752503144732032,0.000723929550450816,0.8773293207500719,
0.0007247633038248635,0.8794172021118196,0.0007256035895739056,0.8815209436139331,0.0007264532597045065,0.8836476979083573,0.0007273152216374242,0.8858047781004107,0.0007281924307543494,0.8879996422658111,
0.0007290878789830281,0.8902398684083558,0.0007300045789273423,0.8925331185794968,0.0007309455430526645,0.8948870908635684,0.0007319137574656844,0.8973094579716888,0.0007329121498886438,0.8998077913001904,
0.0007339435515274697,0.9023894695138696,0.0007350106526783481,0.9050615710285062,0.0007361159521134971,0.9078307502075817,0.000737261700537983,0.910703097667765,0.0001812712131805553,0.0002491384321187703,
0.0254379386800392,0.0003011236389389571,0.09881511572468257,0.0003417143692905089,0.1744135601545568,0.0003741093236630491,0.242998689534953,0.0004004876411801612,0.3031829004247381,0.0004223463640195241,
0.3555733525259414,0.0004407355613198239,0.4012169377011952,0.000456408855964047,0.4411502489205138,0.0004699191284978094,0.476278745476339,0.000481680241910951,0.5073597250977332,0.0004920076781544193,
0.5350167430251173,0.0005011458800384139,0.5597612077529058,0.0005092870408844213,0.5820131769183957,0.00051658427549107,0.6021191036853046,0.0005231610232704372,0.620366279950506,0.0005291178753295886,
0.6369943582199482,0.0005345376082578074,0.65220447386969,0.0005394889485099423,0.6661664586032026,0.0005440294242631437,0.6790245567662841,0.0005482075518868509,0.6909019734144254,0.0005520645308096385,
0.7019045110624699,0.0005556355707323485,0.7121234938350971,0.0005589509407679369,0.7216381322017054,0.0005620368060492128,0.730517446427696,0.0005649159003089793,0.7388218400819079,0.000567608070713959,
0.7466043945017325,0.0005701307223629278,0.7539119395110627,0.000572499183351742,0.7607859437362658,0.0005747270064854298,0.7672632586796091,0.0005768262201098143,0.7733767436165075,0.000578807537812054,
0.7791557928788411,0.0005806805346667902,0.7846267827935451,0.000582453796114398,0.789813452179805,0.0005841350443286614,0.7947372276555162,0.0005857312459743303,0.7994175029023208,0.000587248704505127,
0.8038718793655162,0.0005886931395613307,0.8081163745265718,0.0005900697555567861,0.8121656028098705,0.000591383301170653,0.816032933316133,0.0005926381211586964,0.8197306278697833,0.0005938382016564763,
0.8232699622928324,0.0005949872099502598,0.8266613333475188,0.0005960885295313444,0.8299143534033866,0.0005971452911184121,0.8330379345654622,0.00059816040022474,0.8360403637358284,0.0005991365617580406,
0.8389293698608987,0.0006000763020667812,0.8417121844328502,0.00060098198878518,0.8443955961593316,0.0006018558487773787,0.8469860005853584,0.0006026999844377125,0.849489445340874,0.0006035163885670166,
0.8519116715931933,0.0006043069580133006,0.8542581522025076,0.0006050735062378763,0.8565341270083502,0.0006058177749442889,0.8587446356133339,0.0006065414448864817,0.860894547975822,0.0006072461459539104,
0.862988593073986,0.0006079334666143073,0.865031385858623,0.0006086049627790323,0.8670274526700236,0.0006092621661410419,0.8689812552540616,0.0006099065920211252,0.8708972134736309,0.0006105397467438793,
0.8727797267727406,0.0006111631345506713,0.8746331944112725,0.0006117782640423093,0.8764620344479443,0.000612386654129138,0.8782707014068333,0.0006129898394506246,0.880063702518437,0.000613589375210095,
0.8818456123793366,0.0006141868413530917,0.8836210858249036,0.0006147838459998774,0.8853948687571734,0.0006153820280240385,0.8871718066152805,0.0006159830586501987,0.8889568501193102,0.0006165886419249452,
0.8907550578610847,0.0006172005138967517,0.8925715952587397,0.0006178204403237505,0.8944117293380552,0.0006184502127136527,0.8962808187550784,0.0006190916424892084,0.8981842984351409,0.0006197465530668733,
0.9001276581771762,0.0006204167696375895,0.9021164145644402,0.0006211041064488501,0.9041560755392406,0.0006218103514087189,0.906252097046702,0.0006225372478675618,0.9084098312380428,0.0006232864734841831,
0.9106344658545082,0.0006240596161519326,0.91293095459586,0.000624858147048696,0.9153039385179171,0.0006256833909833134,0.9177576588059725,0.00015202142377128,0.0002089413039628107,0.02544061844919076,
0.0002525434029450244,0.09882898388671699,0.0002865906408583821,0.1744441250170073,0.0003137652564915279,0.2430497235076716,0.0003358945562695622,0.3032570714471141,0.0003542338804453592,0.3556725894600452,
0.0003696637616173895,0.4013426546343276,0.0003828161474542172,0.4413034887132714,0.0003941546527609132,0.4764602748670496,0.0004040263187063507,0.5075700996286872,0.0004126956886076392,0.5352563525449359,
0.0004203677355440847,0.5600303090493647,0.0004272036184847955,0.5823119175989039,0.0004333317268557089,0.6024475399611718,0.0004388555655790971,0.6207243900791432,0.0004438594799484181,0.637382052787274,
0.0004484128767232627,0.652621603723488,0.0004525733807537433,0.6666128210012832,0.0004563892263945745,0.6794999001357331,0.0004599010909485004,0.6914060010021446,0.0004631435158655601,0.7024368836752101,
0.0004661460196359733,0.7126838318080794,0.0004689339774940904,0.7222260166914516,0.0004715293228925426,0.7311324200920039,0.0004739511114194688,0.739463407184761,0.0004762159775818406,0.7472720204610868,
0.0004783385074390325,0.7546050498913561,0.0004803315446135938,0.7615039226740586,0.0004822064431621426,0.768005446719198,0.0004839732777639475,0.7741424349221077,0.0004856410194012131,0.7799442317805824,
0.0004872176829669196,0.78543715961578,0.0004887104519025023,0.7906448982918555,0.0004901257839368306,0.7955888096768238,0.0004914695011954675,0.8002882159858926,0.0004927468673202736,0.8047606394754871,
0.0004939626537433801,0.8090220096175752,0.0004951211968659503,0.8130868428076891,0.0004962264475780244,0.8169683987907728,0.0004972820143036731,0.8206788172836932,0.0004982912005523596,0.8242292376985011,
0.0004992570377925784,0.8276299044001774,0.0005001823143295754,0.8308902595461131,0.0005010696007591025,0.834019025235724,0.000501921272478859,0.8370242764345617,0.0005027395296647642,0.8399135059177911,
0.0005035264150574445,0.8426936822947919,0.0005042838298529326,0.8453713020233193,0.0005050135479486603,0.8479524361928227,0.0005057172287598274,0.8504427727478738,0.0005063964287909282,0.8528476547306911,
0.0005070526121216016,0.8551721150436168,0.0005076871599442239,0.8574209081657413,0.0005083013792721207,0.8595985392007632,0.0005088965109213797,0.8617092905839876,0.0005094737368555486,0.863757246733787,
0.0005100341869706155,0.8657463168957334,0.0005105789453872895,0.8676802563950228,0.0005111090563084296,0.8695626864839279,0.000511625529491321,0.8713971129451628,0.000512129345377131,0.8731869435885925,
0.0005126214599131586,0.8749355047571672,0.000513102809097241,0.8766460569378332,0.0005135743132677963,0.8783218095540483,0.0005140368811573158,0.8799659349980495,0.0005144914137215976,0.8815815819428285,
0.0005149388077515203,0.883171887955587,0.000515379959268618,0.884739991415962,0.0005158157667000665,0.8862890427233201,0.0005162471338228604,0.8878222147576811,0.0005166749724609052,0.8893427125382002,
0.000517100204912442,0.890853782001469,0.000517523766078647,0.892358717799149,0.0005179466052574159,0.8938608699906166,0.0005183696875592966,0.8953636494814934,0.0005187939948953462,0.8968705320333606,
0.0005192205264794681,0.8983850606439769,0.0005196502987807066,0.8999108460714149,0.0005200843448542496,0.9014515652504507,0.0005205237129738048,0.9030109573261521,0.0005209694644829235,0.9045928170091283,
0.0005214226707791797,0.9062009849407291,0.0005218844093433793,0.9078393347463627,0.0005223557587267653,0.9095117564530361,0.0005228377924131475,0.9112221359554993,0.0005233315714807396,0.9129743302365174,
0.0005238381360009637,0.9147721380834732,0.0005243584951293212,0.9166192660984713,0.0005248936158673047,0.918519289874983,0.0005254444105047963,0.920475610313189,0.0005260117227898333,0.922491405170291,
0.0005265963129170898,0.9245695760920818,0.0001268446182127322,0.0001743404286829759,0.02544292517172434,0.0002107251564740404,0.09884092351738187,0.0002391381657362996,0.1744704449940017,0.0002618172499590152,
0.2430936807628649,0.0002802869400172751,0.3033209758285258,0.0002955945748575447,0.355758118174555,0.0003084747892467032,0.4014510452519073,0.0003194548214120799,0.4414356626653468,0.0003289214677172303,
0.4766169192502833,0.0003371642657758759,0.5077517230370061,0.0003444039242024864,0.5354633241883532,0.0003508114513703712,0.5602628868037312,0.0003565213009364801,0.5825702685336475,0.0003616405864628345,
0.6027317556351184,0.0003662556600075083,0.6210344981187766,0.0003704368884625953,0.6377180257702407,0.0003742421752675708,0.652983366071077,0.0003777195940186037,0.667000254371894,0.0003809093836500495,
0.6799128476679434,0.0003838454780905802,0.6918442706570039,0.0003865566919752524,0.7029002508578049,0.000389067649129971,0.7131720413986303,0.0003913995164996082,0.7227387845799725,0.0003935705893726759,
0.7316694342837096,0.0003955967618362022,0.7400243285231667,0.0003974919078428414,0.747856483000171,0.0003992681920658595,0.755212660935947,0.0004009363251645967,0.7621342624977909,0.0004025057747089784,
0.7686580679613385,0.000403984940487517,0.7748168616577529,0.0004053813010180153,0.7806399582528606,0.00040670153662997,0.786153648613477,0.0004079516333750262,0.7913815791511374,0.0004091369711617666,
0.7963450758813084,0.00041026239884156,0.8010634223351132,0.0004113322984474455,0.8055540987877207,0.0004123506403741368,0.8098329889290022,0.000413321030958808,0.813914559025874,0.0004142467536602034,
0.8178120137564664,0.0004151308048232806,0.821537432190936,0.0004159759248469224,0.825101886818934,0.0004167846254347061,0.8285155480533462,0.0004175592134966781,0.8317877762533697,0.0004183018121784048,
0.8349272029911109,0.0004190143794182329,0.8379418030218164,0.0004196987243715124,0.8408389581983645,0.0004203565219890266,0.8436255143875691,0.0004209893259940189,0.8463078322926154,0.0004215985804664407,
0.848891832957265,0.0004221856302130706,0.8513830386190436,0.0004227517300769591,0.8537866094869954,0.0004232980533183947,0.8561073769418921,0.0004238256991815968,0.8583498735907137,0.000424335699746057,
0.8605183605508547,0.0004248290261484384,0.8626168522912827,0.0004253065942498067,0.8646491393164755,0.0004257692698134273,0.8666188089433132,0.000426217873250141,0.8685292643902939,0.0004266531839812305,
0.8703837423717155,0.0004270759444625287,0.8721853293661955,0.0004274868639081453,0.8739369767085355,0.0004278866217474832,0.8756415146360238,0.0004282758708450621,0.8773016654044247,0.0004286552405089854,
0.8789200555747797,0.0004290253393105878,0.8804992275594572,0.0004293867577348235,0.8820416505043871,0.0004297400706782409,0.8835497305738722,0.0004300858398088763,0.8850258206945885,0.0004304246158000579,
0.8864722298061936,0.0004307569404478835,0.887891231657202,0.0004310833486799939,0.8892850731763069,0.0004314043704611773,0.8906559824410115,0.0004317205325992682,0.8920061762571498,0.0004320323604527339,
0.8933378673545269,0.000432340379539238,0.8946532711953979,0.0004326451170423202,0.8959546123837375,0.0004329471032111149,0.8972441306541727,0.000433246872645738,0.8985240864099954,0.0004335449654585917,
0.8997967657697997,0.0004338419282993692,0.901064485071994,0.0004341383152289965,0.9023295947757297,0.0004344346884251275,0.9035944826856999,0.0004347316186991527,0.9048615764168997,0.0004350296858020064,
0.9061333450039083,0.0004353294784934289,0.9074122995477562,0.0004356315943468108,0.9087009927822184,0.0004359366392594052,0.9100020174307473,0.0004362452266356344,0.9113180032156295,0.000436557976209574,
0.9126516123728105,0.0004368755124715989,0.9140055335197583,0.0004371984626638003,0.915382473720423,0.0004375274543093129,0.9167851485915539,0.0004378631122423307,0.9182162702992493,0.000438206055108562,
0.9196785333045526,0.0004385568913104083,0.921174597733176,0.0004389162143774648,0.9227070702680383,0.0004392845977512571,0.9242784824952105,0.0004396625889836182,0.9258912666749617,0.0004400507033609156,
0.9275477289606028,0.0004404494169815118,0.9292500201491721,0.000440859159331359,0.9310001041198086,0.0001055644989267907,0.0001450940149574984,0.02544487491553157,0.0001753772719352769,0.09885101688865597,
0.0001990267012200822,0.1744926990200333,0.0002179045108524649,0.2431308553318966,0.000233279364173155,0.303375032942641,0.0002460228232089687,0.3558304876736915,0.0002567462279756849,0.4015427879008645,
0.0002658883294937232,0.4415475744961243,0.0002737710071803561,0.4767496006654743,0.0002806352044648509,0.5079056259634127,0.0002866645871705139,0.5356387851068741,0.0002920014619138971,0.5604601504270702,
0.0002967577154420555,0.5827895051312842,0.0003010224827883186,0.6029730739429139,0.0003048676217987417,0.6212979551152507,0.0003083516878734,0.6380036341864078,0.0003115228646427098,0.653291100231385,
0.0003144211555866268,0.6673300547852067,0.0003170800443702226,0.6802646246666912,0.0003195277677802438,0.6922179072993687,0.0003217883024392515,0.7032956052494281,0.0003238821374603807,0.7135889485487379,
0.0003258268851960003,0.723177057876695,0.0003276377682377606,0.732128866651248,0.0003293280109066263,0.740504693306939,0.000330909156355051,0.7483575346135667,0.0003323913252386479,0.7557341352926834,
0.0003337834281257755,0.7626758772463315,0.0003350933410056701,0.7692195225318236,0.0003363280511552274,0.7753978371270195,0.0003374937790390407,0.7812401170290845,0.0003385960807104609,0.7867726339385476,
0.0003396399342555158,0.792019014415911,0.0003406298131057851,0.7970005637462622,0.0003415697484891286,0.8017365436465392,0.0003424633828504657,0.8062444112774049,0.0003433140157303645,0.8105400256832551,
0.0003441246433158771,0.8146378267077895,0.0003448979926599179,0.8185509905633486,0.0003456365513904324,0.82229156552692,0.0003463425935893944,0.8258705906609294,0.0003470182024071915,0.8292981999865211,
0.0003476652898847084,0.8325837141504723,0.0003482856143791181,0.835735721307981,0.000348880795926679,0.8387621486794669,0.0003494523298240878,0.8416703260200142,0.0003500015986670602,0.8444670420569842,
0.000350529883049145,0.8471585947980717,0.0003510383710940082,0.8497508364833812,0.0003515281669694799,0.8522492138466825,0.0003520002985106894,0.8546588042593751,0.000352455724061933,0.8569843482530394,
0.0003528953386319557,0.859230278850431,0.0003533199794446321,0.8614007480784993,0.0003537304309562248,0.8634996509889044,0.0003541274294011706,0.8655306474702691,0.000354511666920449,0.8674971821009717,
0.0003548837953198046,0.8694025022607405,0.0003552444294992589,0.8712496746929433,0.0003555941505903093,0.8730416006866109,0.0003559335088328495,0.8747810300273897,0.0003562630262200606,0.8764705738493466,
0.0003565831989362264,0.8781127165044531,0.0003568944996095447,0.8797098265533606,0.0003571973793994888,0.8812641669694672,0.0003574922699360495,0.8827779046380256,0.0003577795851262342,0.8842531192229806,
0.0003580597228414603,0.8856918114661521,0.0003583330664979325,0.8870959109761732,0.0003585999865407079,0.888467283558108,0.0003588608418408972,0.8898077381288096,0.0003591159810143158,0.891119033257723,
0.0003593657436688535,0.8924028833679229,0.0003596104615868662,0.8936609646275998,0.000359850459847991,0.8948949205579224,0.0003600860578969333,0.8961063673791346,0.0003603175705599569,0.8972968991128309,
0.0003605453090130204,0.8984680924545586,0.0003607695817037281,0.8996215114271492,0.0003609906952284962,0.9007587118214661,0.0003612089551655714,0.9018812454275019,0.0003614246668637566,0.9029906640549621,
0.0003616381361859153,0.9040885233385761,0.0003618496702055102,0.9051763863193673,0.0003620595778536006,0.9062558267889668,0.0003622681705128622,0.9073284323797515,0.0003624757625543019,0.9083958073791226,
0.0003626826718114261,0.90945957524161,0.0003628892199856792,0.9105213807677069,0.0003630957329760065,0.9115828919134202,0.0003633025411244249,0.9126458011895109,0.0003635099793685103,0.9137118266043415,
0.0003637183872907533,0.9147827120992092,0.0003639281090538103,0.9158602274201335,0.0003641394932098118,0.9169461673653735,0.0003643528923711136,0.9180423503436506,0.0003645686627292243,0.9191506161742795,
0.0003647871634081525,0.9202728230574015,0.0003650087556381453,0.9214108436404823,0.0003652338017357785,0.9225665601064491,0.0003654626638766759,0.9237418582096152,0.0003656957026478399,0.9249386201881679,
0.0003659332753677379,0.9261587164868569,0.0003661757341639742,0.9274039962309573,0.000366423423800659,0.9286762764029628,0.0003666766792505202,0.9299773296871511,0.0003669358230104597,0.9313088709644671,
0.0003672011621636608,0.9326725424613703,0.0003674729851965463,0.9340698975815723,0.000367751558584863,0.9355023834790327,0.0003680371231699034,0.9369713224641266,8.780076664716277e-05,0.0001206797894237555,
0.02544650247967021,0.0001458689907141319,0.09885944344034912,0.0001655410204684315,0.1745112808102761,0.0001812446142239127,0.243161901103803,0.0001940348805295792,0.3034201872483507,0.0002046366722127249,
0.3558909526523158,0.0002135584233839305,0.4016194593852736,0.0002211650307991108,0.4416411290854953,0.0002277241922738999,0.4768605533838495,0.0002334362932843276,0.5080343695969284,0.0002384540829843639,
0.5357856177524436,0.0002428959131113313,0.5606252952358727,0.0002468548359474569,0.5829731247513873,0.0002504049817295676,0.6031752812789858,0.0002536061116756343,0.6215188214764482,0.0002565069236788439,
0.63824319554714,0.0002591474896792285,0.6535493621049741,0.0002615610783797454,0.6676069960615446,0.0002637755361065909,0.6805602006643017,0.0002658143454788651,0.692532052212964,0.0002676974460335837,
0.70362823412551,0.0002694418768222108,0.7139399588872692,0.0002710622843532818,0.723546330931951,0.0002725713276157303,0.732516268486713,0.0002739800036682437,0.740910075645583,0.0002752979113613955,
0.7487807355146596,0.0002765334664638627,0.7561749796785543,0.0002776940783128682,0.7631341772964979,0.0002787862957738123,0.7696950779572316,0.0002798159285470986,0.7758904353334187,0.0002807881485415345,
0.7817495331755639,0.0002817075750299807,0.7872986308948055,0.0002825783465328497,0.7925613426197803,0.000283404181779784,0.7975589609612814,0.0002841884316364275,0.8023107346178519,0.0002849341235200138,
0.8068341072829448,0.0002856439995410211,0.8111449239759836,0.0002863205493808272,0.815257609843654,0.0002869660387338721,0.8191853256086143,0.0002875825339972447,0.8229401031375546,0.0002881719237731651,
0.8265329640257898,0.0002887359376546216,0.8299740236251799,0.0002892761626868619,0.8332725825556256,0.000289794057833975,0.8364372074214862,0.0002902909667276428,0.8394758021901656,0.0002907681289320944,
0.8423956714705967,0.0002912266899236334,0.8452035767462418,0.0002916677099534407,0.847905786463963,0.0002920921719375913,0.8505081207514051,0.0002925009884974771,0.8530159914271018,0.0002928950082563857,
0.8554344378758767,0.0002932750214832763,0.8577681592844462,0.0002936417651623478,0.8600215436661002,0.0002939959275564363,0.8621986940470617,0.0002943381523232877,0.8643034521390155,0.0002946690422360806,
0.866339419781072,0.0002949891625530086,0.8683099783990099,0.0002952990440750963,0.8702183066991327,0.0002955991859265754,0.8720673967877357,0.0002958900580879658,0.8738600688843806,0.0002961721037083923,
0.8755989847774061,0.0002964457412205325,0.8772866601529093,0.0002967113662788715,0.8789254759134587,0.0002969693535395673,0.880517688589716,0.0002972200582981628,0.8820654399367039,0.0002974638179995708,
0.8835707657964145,0.0002977009536331723,0.8850356042996281,0.0002979317710244726,0.886461803472041,0.0002981565620335326,0.8878511283029286,0.0002983756056693069,0.8892052673284969,0.0002985891691280583,
0.8905258387766696,0.0002987975087631656,0.8918143963152597,0.0002990008709928799,0.8930724344411763,0.0002991994931519033,0.8943013935444719,0.0002993936042920544,0.8955026646775683,0.0002995834259367309,
0.896677594056873,0.0002997691727933844,0.8978274873211489,0.0002999510534277608,0.8989536135684026,0.0003001292709032479,0.9000572091906636,0.0003003040233882832,0.901139481523815,0.0003004755047344127,
0.9022016123275653,0.0003006439050272616,0.9032447611087047,0.0003008094111123485,0.9042700682989381,0.0003009722070973778,0.9052786582968102,0.0003011324748323426,0.9062716423815145,0.0003012903943684867,
0.9072501215046918,0.000301446144396887,0.9082151889646557,0.0003015999026671373,0.9091679329658161,0.0003017518463863298,0.9101094390643967,0.0003019021525982453,0.9110407924998413,0.00030205099854237,
0.9119630804095697,0.0003021985619920632,0.9128773939229607,0.0003023450215708908,0.9137848301286094,0.0003024905570458278,0.9146864939070101,0.0003026353495957124,0.9155834996188661,0.0003027795820529983,
0.9164769726372038,0.0003029234391165164,0.917368050709395,0.0003030671075326073,0.9182578851330556,0.0003032107762416357,0.9191476417276118,0.0003033546364865421,0.9200385015811228,0.0003034988818797365,
0.9209316615497335,0.0003036437084242899,0.9218283344849421,0.0003037893144850508,0.9227297491617354,0.0003039359007049947,0.923637149878606,0.0003040836698618365,0.9245517956985898,0.0003042328266596891,
0.9254749592987896,0.000304383577450359,0.9264079253944743,0.0003045361298787459,0.9273519887028263,0.0003046906924467665,0.9283084514108677,0.0003048474739902782,0.9292786201121033,0.0003050066830636484,
0.9302638021771222,0.0003051685272269233,0.931265301524892,0.0003053332122310158,0.9322844137639074,0.0003055009410969733,0.9333224206758433,0.0003056719130862303,0.9343805840190375,0.0003058463225598092,
0.9354601386351341,0.0003060243577257252,0.9365622848496463,0.0003062061992753964,0.9376881801661726,0.0003063920189116574,0.9388389302645763,0.0003065819777730374,0.940015579325664,0.0003067762247612739,
0.941219099718771,0.0003069748947815949,0.9424503811041136,7.309194121921338e-05,0.0001004637913909749,0.02544785012683798,0.0001214344266344417,0.09886642145889683,0.0001378124012607418,0.174526670295946,
0.0001508869219533374,0.2431876171096405,0.0001615362906618204,0.3034575961449956,0.0001703638961634086,0.3559410557392941,0.0001777929742443069,0.4016830055695095,0.0001841272759437895,0.4417186869611078,
0.0001895896409160979,0.4769525589580654,0.0001943468759030283,0.5081411588464865,0.0001985261336437941,0.5359074494195392,0.0002022259333499051,0.5607623671619424,0.0002055237341359384,0.5831255858225246,
0.0002084812437801258,0.6033432403229687,0.0002111482088116174,0.6217023539986306,0.0002135651662786688,0.6384423489057279,0.0002157654726972967,0.6537641595402491,0.0002177768213406241,0.6678374398703552,
0.0002196223917129092,0.6808062747278515,0.0002213217308213229,0.692793724028602,0.0002228914362913128,0.7039054564531154,0.0002243456912849616,0.7142326710923299,0.0002256966873287446,0.7238544600868844,
0.0002269549614676161,0.7328397302769779,0.0002281296672953758,0.7412487751160227,0.0002292287944844476,0.7491345676832047,0.0002302593478625398,0.7565438300381137,0.0002312274944605136,0.7635179222210341,
0.0002321386850119108,0.7700935850242779,0.000232997754930363,0.7763035635723555,0.0002338090086934327,0.7821771332486441,0.0002345762907259228,0.7877405452160454,0.0002353030452346481,0.7930174054153138,
0.0002359923669511411,0.7980289982735111,0.0002366470443530013,0.8027945642546851,0.0002372695966322588,0.807331538712508,0.000237862305440656,0.8116557581664309,0.0002384272422525216,0.8157816390470166,
0.0002389662920348844,0.8197223330870221,0.0002394811737932779,0.823489862829619,0.0002399734584639214,0.8270952401504242,0.0002404445845437014,0.8305485702196599,0.0002408958717848132,0.8338591429442381,
0.0002413285332280931,0.8370355136106796,0.0002417436858056495,0.8400855741856877,0.0002421423597075709,0.8430166165116897,0.0002425255066777939,0.845835388451546,0.0002428940073795231,0.8485481438833696,
0.0002432486779499734,0.8511606873176778,0.0002435902758469355,0.8536784138006696,0.0002439195050751448,0.8561063446757747,0.0002442370208681859,0.8584491596979498,0.0002445434338913074,0.8607112259291628,
0.0002448393140217253,0.8628966237872217,0.0002451251937555122,0.8650091705719897,0.0002454015712837818,0.8670524417517986,0.0002456689132754098,0.8690297902574413,0.0002459276573988435,0.8709443640006193,
0.0002461782146115141,0.8727991218073733,0.0002464209712418874,0.8745968479342314,0.0002466562908861773,0.8763401653150386,0.0002468845161391411,0.87803154766925,0.0002471059701761088,0.8796733305874987,
0.0002473209582014312,0.8812677216971869,0.0002475297687768076,0.8828168099994235,0.0002477326750414573,0.884322574458617,0.0002479299358347776,0.8857868919172412,0.0002481217967309792,0.8872115444005592,
0.0002483084909941748,0.8885982258692703,0.0002484902404615005,0.8899485484720331,0.0002486672563610595,0.8912640483444852,0.0002488397400707827,0.892546190996671,0.0002490078838236775,0.8937963763265875,
0.0002491718713643894,0.8950159432938443,0.0002493318785615123,0.896206174284098,0.0002494880739796445,0.8973682991919631,0.0002496406194148037,0.8985034992474406,0.0002497896703964584,0.8996129106085194,
0.0002499353766591264,0.9006976277404708,0.0002500778825862051,0.9017587066004092,0.000250217327628446,0.9027971676439548,0.0002503538466992535,0.9038139986692414,0.0002504875705487804,0.9048101575120641,
0.0002506186261185988,0.9057865746046461,0.0002507471368785526,0.90674415540928,0.000250873223147232,0.907683782736985,0.000250997002397365,0.9086063189602724,0.0002511185895472782,0.9095126081281425,
0.000251238097239448,0.9104034779905169,0.0002513556361070412,0.9112797419384429,0.0002514713150292236,0.9121422008655768,0.0002515852413759041,0.9129916449556555,0.0002516975212424692,0.9138288553998849,
0.0002518082596749575,0.9146546060474243,0.0002519175608860183,0.915469664991386,0.0002520255284618921,0.9162747960920304,0.0002521322655605479,0.917070760438088,0.0002522378751010061,0.9178583177463857,
0.0002523424599437699,0.9186382276991902,0.000252446123062179,0.9194112512179023,0.0002525489677043914,0.9201781516709387,0.0002526510975455835,0.9209396960128159,0.0002527526168298444,0.9216966558506134,
0.0002528536305011233,0.922449808433123,0.0002529542443224635,0.9231999375571039,0.0002530545649826379,0.923947834384147,0.0002531547001891691,0.9246942981607159,0.0002532547587465925,0.9254401368329844,
0.0002533548506186903,0.9261861675471152,0.0002534550869732947,0.9269332170246625,0.0002535555802081323,0.9276821218018076,0.0002536564439560544,0.928433728320188,0.0002537577930678808,0.9291888928561548,
0.0002538597435709706,0.9299484812744166,0.0002539624126015323,0.93071336859122,0.0002540659183085976,0.9314844383314878,0.0002541703797275094,0.9322625816637383,0.0002542759166207297,0.9330486962961403,
0.0002543826492837433,0.9338436851167795,0.000254490698313846,0.9346484545611476,0.0002546001843396463,0.9354639126900528,0.0002547112277091971,0.9362909669616406,0.0002548239481348102,0.9371305216820521,
0.0002549384642927942,0.9379834751204702,0.0002550548933766077,0.9388507162759753,0.0002551733506022396,0.9397331212857866,0.0002552939486650164,0.9406315494671578,0.0002554167971475084,0.941546838988458,
0.0002555420018787558,0.9424798021688484,0.0002556696642456735,0.9434312204104797,0.0002557998804582168,0.9444018387723074,0.0002559327407706962,0.9453923602004483,0.0002560683286625224,0.9464034394364762,
0.0002562067199826295,0.9474356766321346,6.097078489116776e-05,8.380405076108355e-05,0.0254489606593285,0.0001012978969826672,0.09887217219781226,0.0001149608733516808,0.1745393544210289,0.0001258683390053939,
0.2432088150280462,0.0001347529092969161,0.3034884369920806,0.0001421178877684437,0.3559823686685004,0.0001483163075228666,0.4017354125928528,0.0001536015360676195,0.4417826623358451,0.0001581594503691672,
0.4770284685054013,0.0001621291886747866,0.5082292867703224,0.0001656168126296357,0.5360080170824175,0.0001687044997776278,0.5608755465361002,0.0001714568607858636,0.5832515093250096,0.0001739253676423436,
0.6034820081391302,0.0001761515150836057,0.6218540396252812,0.0001781691159294222,0.6386070034090316,0.0001800059934918693,0.6539418148701861,0.0001816852471937282,0.6680281114720945,0.0001832262113807821,
0.6810099636217561,0.0001846451904182646,0.6930104184833876,0.0001859560284988198,0.7041351333436702,0.0001871705558346453,0.7144752970102477,0.0001882989413512779,0.7241099922548785,0.0001893499739181992,
0.733108117306014,0.0001903312884234792,0.7415299576349031,0.000191249548890067,0.7494284788635305,0.0001921105978487627,0.7568503960324043,0.0001929195789948627,0.7638370625278608,0.0001936810385340217,
0.7704252127912312,0.0001943990094098548,0.7766475858453091,0.0001950770816902051,0.7825334511739168,0.000195718461692066,0.7881090542005847,0.0001963260218904405,0.7933979952488399,0.0001969023432430833,
0.7984215532155936,0.0001974497512412977,0.8031989630889407,0.0001979703467447663,0.8077476547694435,0.000198466032459483,0.8120834593159135,0.0001989385357600094,0.8162207876608769,0.0001993894284313026,
0.8201727859718991,0.0001998201438042677,0.8239514711297961,0.0002002319916776369,0.8275678492200949,0.0002006261713526632,0.831032019463781,0.0002010037830532616,0.8343532656268658,0.0002013658379601603,
0.8375401366294501,0.000201713267051409,0.8406005178108765,0.0002020469289116992,0.8435416940880704,0.0002023676166481832,0.8463704060610679,0.0002026760640298827,0.849092899966467,0.0002029729509505747,
0.851714972250835,0.0002032589083006388,0.8542420094276686,0.0002035345223212339,0.8566790237898599,0.0002038003385039595,0.859030685471952,0.0002040568650905127,0.8613013512904316,0.0002043045762195164,
0.8634950907340174,0.0002045439147614532,0.865615709427792,0.000204775294877309,0.8676667703537839,0.0002049991043319707,0.8696516130751836,0.0002052157065895086,0.8715733711808579,0.0002054254427141078,
0.8734349881404843,0.000205628633097507,0.8752392317378276,0.0002058255790312958,0.876988707229908,0.0002060165641402434,0.8786858693626248,0.0002062018556909442,0.8803330333584298,0.0002063817057884218,
0.8819323849785773,0.0002065563524718971,0.8834859897510554,0.0002067260207196752,0.8849958014452854,0.0002068909233720063,0.8864636698658906,0.0002070512619798147,0.8878913480301039,0.0002072072275863421,
0.8892804987865693,0.0002073590014480071,0.8906327009272852,0.0002075067557001234,0.8919494548391176,0.0002076506539725411,0.8932321877366053,0.0002077908519597573,0.8944822585136047,0.0002079274979495885,
0.8957009622476075,0.0002080607333140918,0.8968895343872672,0.0002081906929660589,0.8980491546507244,0.0002083175057840894,0.8991809506596934,0.00020844129500896,0.9002860013319247,0.0002085621786137507,
0.9013653400525538,0.00020868026964996,0.9024199576429591,0.0002087956765716353,0.9034508051440568,0.0002089085035393573,0.9044587964294366,0.0002090188507057549,0.9054448106623624,0.000209126814484072,
0.906409694609428,0.000209232487801176,0.9073542648225251,0.0002093359603362731,0.9082793096997696,0.0002094373187464826,0.9091855914351026,0.0002095366468803258,0.910073847865441,0.0002096340259800857,
0.9109447942234818,0.0002097295348739169,0.9117991248035615,0.0002098232501585007,0.9126375145473211,0.0002099152463729734,0.9134606205553339,0.0002100055961647897,0.9142690835302997,0.000210094370448118,
0.9150635291568929,0.0002101816385553131,0.9158445694228796,0.0002102674683819544,0.9166128038856619,0.0002103519265258884,0.9173688208879894,0.0002104350784206675,0.9181131987261704,0.0002105169884637311,
0.9188465067737359,0.0002105977201396311,0.919569306563135,0.0002106773361385627,0.9202821528276868,0.0002107558984704203,0.9209855945056615,0.0002108334685745573,0.9216801757080246,0.0002109101074253908,
0.9223664366510376,0.0002109858756339512,0.9230449145545729,0.0002110608335454401,0.9237161445066656,0.0002111350413328171,0.9243806602944845,0.0002112085590864012,0.9250389952015636,0.0002112814468994273,
0.9256916827707871,0.0002113537649494626,0.9263392575322679,0.0002114255735755424,0.9269822556948962,0.0002114969333508442,0.9276212157999659,0.0002115679051506757,0.9282566793349052,0.0002116385502155082,
0.9288891913047499,0.0002117089302087404,0.9295193007585985,0.0002117791072688329,0.9301475612678799,0.0002118491440554076,0.9307745313528513,0.0002119191037888559,0.9314007748533154,0.0002119890502829553,
0.9320268612391253,0.0002120590479699442,0.9326533658556098,0.0002121291619174575,0.9332808700986251,0.0002121994578366795,0.9339099615135131,0.0002122700020810254,0.9345412338118337,0.0002123408616346198,
0.9351752867993358,0.0002124121040897991,0.9358127262082568,0.0002124837976128301,0.9364541634266929,0.000212556010897005,0.9371002151174705,0.0002126288131022461,0.9377515027186955,0.0002127022737803375,
0.9384086518179469,0.0002127764627848895,0.9390722913919686,0.0002128514501651432,0.9397430529036528,0.0002129273060427348,0.9404215692481872,0.0002130041004705666,0.9411084735403966,0.0002130819032729711,
0.9418043977356299,0.0002131607838664165,0.9425099710769944,0.0002132408110600791,0.9432258183623697,0.0002133220528357081,0.9439525580254419,0.0002134045761063316,0.9446908000260154,0.0002134884464535012,
0.9454411435460905,0.000213573727842945,0.9462041744896647,0.0002136604823187037,0.9469804627859271,0.0002137487696760521,0.9477705594974885,0.0002138386471137701,0.9485749937375257,0.0002139301688666145,
0.9493942694022146,0.0002140233858191606,0.9502288617275928,0.0002141183451025241,0.9510792136830037,0.0002142150896758442,0.9519457322165201,5.038843819242724e-05,6.925910573206593e-05,0.02544993017457056,
8.3717269029937e-05,0.09887719305509406,9.500958703765662e-05,0.1745504296500425,0.0001040247370141086,0.2432273260547244,0.0001113681540690559,0.3035153719186541,0.0001174557699843303,0.3560184543414467,
0.0001225793327188225,0.4017811956457096,0.0001269482262440711,0.4418385610275083,0.000130716060065286,0.477094807039437,0.0001339978185876816,0.5083063187019734,0.0001368811515412669,0.5360959417043054,
0.0001394339707035654,0.5609745206153379,0.0001417096704842709,0.5833616556836523,0.0001437507874107348,0.6036034220420834,0.0001455916127490372,0.6219867936141575,0.0001472600893816645,0.6387511510142359,
0.0001487792104243456,0.6540973934997754,0.0001501680651407632,0.6681951446873501,0.0001514426313109642,0.6811884629495294,0.0001526163827196775,0.6932003848740875,0.0001537007600492208,0.7043365583549521,
0.0001547055396155196,0.7146881637760102,0.0001556391248364784,0.7243342762844464,0.0001565087786427017,0.7333437871476419,0.0001573208103069447,0.7417769754298066,0.0001580807267724731,0.7496868008107239,
0.0001587933560956882,0.7571199727801611,0.0001594629488101686,0.7641178395042019,0.0001600932616793006,0.7707171304832252,0.0001606876273022243,0.7769505800349328,0.0001612490122811637,0.782847453136627,
0.0001617800660822615,0.788433990871492,0.0001622831622801481,0.7937337893603658,0.0001627604335348945,0.7987681234096885,0.0001632138013840797,0.803556224006277,0.000163645001724293,0.808115517117466,
0.0001640556066920045,0.8124618299171829,0.0001644470435232973,0.8166095694827731,0.0001648206108678467,0.8205718781384362,0.0001651774929489843,0.8243607689160389,0.0001655187718942932,0.8279872440294369,
0.0001658454385065429,0.8314613987881412,0.0001661584017002631,0.8347925129896938,0.0001664584967928401,0.8379891315112622,0.0001667464928090865,0.8410591355569131,0.0001670230989335329,0.8440098057975293,
0.000167288970224224,0.8468478784572497,0.000167544712684776,0.8495795952470572,0.0001677908877772391,0.8522107479174399,0.000168028016446403,0.8547467180936171,0.0001682565827161736,0.8571925129651889,
0.0001684770369102029,0.8595527973243928,0.0001686897985418176,0.8618319223811186,0.000168895258912224,0.8640339517265517,0.0001690937834508111,0.8661626847691907,0.0001692857138269685,0.8682216779257586,
0.0001694713698590671,0.8702142638140924,0.0001696510512430158,0.8721435686645843,0.0001698250391200247,0.8740125281403963,0.0001699935975008061,0.8758239017338748,0.0001701569745613678,0.8775802858868165,
0.0001703154038237586,0.8792841259650483,0.0001704691052335605,0.8809377272028146,0.000170618286144567,0.882543264719395,0.0001707631422198987,0.8841027926989548,0.0001709038582577726,0.8856182528146059,
0.0001710406089492363,0.8870914819688771,0.0001711735595743775,0.8885242194150492,0.0001713028666428259,0.889918113317008,0.000171428678483743,0.8912747267992474,0.0001715511357899551,0.892595543533345,
0.0001716703721204033,0.8938819729025228,0.0001717865143646592,0.8951353547817279,0.000171899683172879,0.8963569639669645,0.0001720099933542312,0.8975480142843038,0.0001721175542465397,0.8987096624060582,
0.0001722224700596145,0.8998430113989865,0.000172324840194509,0.9009491140270437,0.0001724247595407294,0.9020289758291009,0.0001725223187532327,0.9030835579901713,0.0001726176045108827,0.9041137800229966,
0.0001727106997578779,0.9051205222753304,0.0001728016839295312,0.9061046282768895,0.0001728906331636581,0.9070669069387198,0.0001729776204987175,0.9080081346166083,0.0001730627160597527,0.9089290570491781,
0.0001731459872330869,0.9098303911803945,0.0001732274988306475,0.9107128268753881,0.0001733073132447196,0.9115770285377645,0.0001733854905938624,0.9124236366358845,0.0001734620888606598,0.9132532691449957,
0.0001735371640219225,0.9140665229115263,0.0001736107701719091,0.9148639749453517,0.0001736829596390858,0.9156461836453693,0.0001737537830969034,0.9164136899633,0.0001738232896690321,0.9171670185102322,
0.0001738915270294578,0.9179066786100761,0.0001739585414978131,0.9186331653037572,0.0001740243781302842,0.9193469603076764,0.0001740890808064089,0.9200485329296795,0.0001741526923120535,0.9207383409455177,
0.0001742152544188343,0.9214168314385361,0.0001742768079602262,0.9220844416050994,0.0001743373929045788,0.922741599528049,0.0001743970484252448,0.9233887249202853,0.0001744558129680012,0.9240262298403793,
0.000174513724315931,0.9246545193819353,0.0001745708196519153,0.9252739923382532,0.000174627135618867,0.9258850418436779,0.0001746827083778257,0.9264880559928549,0.0001747375736640159,0.9270834184389665,
0.0001747917668409568,0.9276715089718623,0.0001748453229526972,0.9282527040768545,0.0001748982767742368,0.9288273774748004,0.0001749506628601784,0.9293959006439499,0.0001750025155916436,0.9299586433238869,
0.0001750538692214704,0.9305159740017549,0.0001751047579176975,0.9310682603808031,0.0001751552158053241,0.9316158698311472,0.0001752052770063228,0.9321591698224861,0.0001752549756778664,0.932698528338364,
0.0001753043460487161,0.9332343142714111,0.0001753534224537021,0.9337668977988398,0.0001754022393662154,0.9342966507373076,0.0001754508314286094,0.9348239468760963,0.0001754992334803999,0.9353491622873846,
0.0001755474805841303,0.9358726756122197,0.0001755956080487566,0.9363948683206198,0.0001756436514503888,0.9369161249440563,0.0001756916466502072,0.937436833278389,0.0001757396298093589,0.9379573845551401,
0.0001757876374006189,0.9384781735788145,0.0001758357062165868,0.9389995988277885,0.0001758838733741719,0.9395220625161108,0.0001759321763151046,0.9400459706133824,0.0001759806528021968,0.9405717328197126,
0.0001760293409110582,0.9410997624925852,0.0001760782790169647,0.9416304765223171,0.0001761275057765599,0.942164295152655,0.0001761770601040649,0.9427016417429337,0.0001762269811416588,0.9432429424681213,
0.0001762773082236904,0.9437886259529981,0.0001763280808343769,0.944339122836676,0.0001763793385586478,0.9448948652636469,0.0001764311210257938,0.9454562862975816,0.0001764834678455926,0.9460238192541732,
0.0001765364185365951,0.9465978969494375,0.0001765900124462755,0.9471789508600694,0.000176644288662773,0.9477674101926942,0.0001766992859179855,0.9483637008591626,0.0001767550424818139,0.9489682443554301,
0.0001768115960474016,0.949581456542025,0.0001768689836072686,0.9502037463246646,0.0001769272413203024,0.9508355142342285,0.0001769864043696395,0.951477150906035,0.0001770465068115518,0.9521290354592113,
0.0001771075814155438,0.9527915337778898,0.0001771696594959625,0.9534649966970047,0.0001772327707355348,0.9541497580966158,0.0001772969430013587,0.954846132909919,0.0001773622021540052,0.9555544150514501,
0.0001774285718505145,0.9562748752734012,4.141439367791228e-05,5.692455280094251e-05,0.02545075230958515,6.880819048293672e-05,0.0988814509305572,7.808989690617372e-05,0.174559822567384,8.550005185877479e-05,
0.2432430266865816,9.153624710259941e-05,0.3035382198764821,9.654032678652593e-05,0.3560490681596453,0.000100752071555041,0.401820041498867,0.0001043435691172141,0.4418859966836131,0.0001074410674826558,
0.4771511109577878,0.0001101390695970056,0.5083717098779069,0.0001125096126676353,0.5361705935601133,0.0001146085117292486,0.5610585709493248,0.0001164796503201761,0.5834552139424402,0.000118157988128934,
0.6037065747337148,0.0001196717082477531,0.6220996084034542,0.0001210437761733641,0.6388736798670074,0.0001222930892987222,0.6542296751290908,0.000123435336526936,0.6683372064799513,0.0001244836495012513,
0.6813403224997011,0.000125449101886255,0.6933620512173756,0.0001263410963843689,0.7045080329697024,0.0001271676677918698,0.7148694414033631,0.0001279357225503386,0.7245253456049399,0.000128651229759923,
0.7335446313434784,0.00012931937473035,0.7419875726551162,0.0001299446833543593,0.7499071245871521,0.0001305311235624521,0.7573499923313203,0.0001310821886317202,0.7643575200396533,0.0001316009660202284,
0.7709664334404412,0.0001320901945745293,0.7772094632859758,0.0001325523123360135,0.7831158711649373,0.0001329894966984371,0.7887118949231211,0.0001334036983057875,0.7940211275731494,0.0001337966697989137,
0.7990648409231811,0.0001341699903007926,0.8038622630547293,0.0001345250863590139,0.8084308171076828,0.0001348632499289629,0.8127863274927396,0.0001351856538739728,0.8169431985757714,0.0001354933653731529,
0.8209145700097322,0.0001357873575589372,0.8247124521846697,0.000136068519651005,0.8283478446917955,0.0001363376658083228,0.8318308402273038,0.000136595542884475,0.8351707159751699,0.0001368428372415219,
0.8383760141893339,0.0001370801807530195,0.8414546134316254,0.0001373081561065382,0.8444137917023125,0.0001375273014991923,0.8472602825170737,0.0001377381148057014,0.8500003248309528,0.0001379410572868236,
0.8526397075811507,0.0001381365568962136,0.855183809512096,0.0001383250112355325,0.8576376348545903,0.0001385067902006957,0.860005845353164,0.0001386822383562781,0.8622927890697443,0.0001388516770701086,
0.8645025263354547,0.0001390154064358504,0.8666388531742509,0.0001391737070077414,0.8687053224808631,0.0001393268413685703,0.8707052632000881,0.0001394750555493091,0.8726417977239584,0.0001396185803165322,
0.8745178576969695,0.0001397576323417808,0.8763361983967451,0.0001398924152653265,0.8780994118377522,0.0001400231206653106,0.8798099387284803,0.0001401499289419495,0.8814700793975387,0.0001402730101253844,
0.8830820037910464,0.0001403925246147742,0.8846477606322743,0.0001405086238553824,0.8861692858244719,0.0001406214509596627,0.8876484101690278,0.0001407311412776921,0.8890868664633782,0.0001408378229217274,
0.8904862960362619,0.0001409416172491532,0.8918482547719127,0.000141042639307643,0.8931742186694573,0.000141140998245961,0.8944655889790835,0.0001412367976934811,0.8957236969523653,0.0001413301361111912,
0.8969498082404223,0.0001414211071166758,0.8981451269702901,0.0001415097997853232,0.8993107995269396,0.0001415962989297887,0.9004479180657531,0.0001416806853595467,0.9015575237779277,0.0001417630361221969,
0.9026406099291728,0.0001418434247280258,0.9036981246901897,0.0001419219213591954,0.9047309737757363,0.0001419985930647974,0.9057400229075616,0.0001420735039429043,0.9067261001151347,0.0001421467153106468,
0.9076899978868622,0.0001422182858632548,0.908632475183385,0.0001422882718229183,0.9095542593235423,0.0001423567270782516,0.9104560477526908,0.0001424237033150753,0.9113385097022479,0.0001424892501391719,
0.9122022877485889,0.0001425534151916154,0.9130479992787559,0.0001426162442572255,0.9138762378698273,0.0001426777813666525,0.9146875745882436,0.0001427380688925579,0.9154825592148789,0.0001427971476403172,
0.916261721401193,0.0001428550569336404,0.917025571761375,0.0001429118346954711,0.9177746029050118,0.0001429675175245004,0.9185092904144644,0.0001430221407676021,0.9192300937708135,0.0001430757385884771,
0.9199374572319491,0.000143128344032768,0.9206318106661012,0.000143179989089888,0.9213135703438758,0.0001432307047517914,0.9219831396916189,0.0001432805210688911,0.9226409100087367,0.0001433294672033191,
0.923287261151399,0.0001433775714797077,0.9239225621848796,0.000143424861433657,0.9245471720066236,0.0001434713638580431,0.9251614399419787,0.000143517104847309,0.9257657063143877,0.0001435621098398706,
0.9263603029917071,0.000143606403658759,0.9269455539101992,0.0001436500105506133,0.9275217755776247,0.000143692954223128,0.9280892775567625,0.0001437352578810521,0.9286483629305795,0.0001437769442608283,
0.9291993287501823,0.0001438180356639552,0.9297424664665896,0.0001438585539891476,0.9302780623472868,0.0001438985207633644,0.9308063978784351,0.0001439379571717663,0.9313277501535416,0.0001439768840866617,
0.9318423922493136,0.0001440153220954916,0.9323505935893585,0.0001440532915279006,0.932852620296316,0.0001440908124819338,0.9333487355329478,0.0001441279048493981,0.933839199832643,0.0001441645883404158,
0.934324271419735,0.0001442008825071996,0.9348042065199632,0.0001442368067670687,0.9352792596613518,0.000144272380424723,0.935749683965716,0.0001443076226937871,0.9362157314309451,0.0001443425527176319,
0.936677653204151,0.0001443771895894746,0.9371356998457099,0.0001444115523717566,0.9375901215841609,0.0001444456601147895,0.9380411685618638,0.0001444795318746588,0.9384890910712563,0.0001445131867303672,
0.9389341397814829,0.0001445466438001951,0.9393765659551052,0.0001445799222572509,0.9398166216545342,0.000144613041344178,0.9402545599377586,0.0001446460203869812,0.9406906350428731,0.000144678878807928,
0.9411251025608388,0.0001447116361374763,0.9415582195958393,0.0001447443120251743,0.9419902449125173,0.0001447769262494714,0.9424214390693088,0.000144809498726375,0.942852064537012,0.0001448420495168817,
0.943282385801656,0.0001448745988331048,0.9437126694506554,0.0001449071670430162,0.9441431842411634,0.0001449397746737126,0.9445742011494633,0.0001449724424131142,0.9450059934001559,0.0001450051911099938,
0.9454388364738392,0.0001450380417722341,0.9458730080918976,0.0001450710155632024,0.946308788176958,0.0001451041337961289,0.9467464587875015,0.0001451374179263711,0.9471863040250684,0.000145170889541443,
0.9476286099124387,0.0001452045703486834,0.9480736642411287,0.0001452384821604381,0.9485217563865096,0.0001452726468766264,0.94897317708883,0.0001453070864645653,0.9494282181984106,0.0001453418229359214,
0.9498871723832844,0.0001453768783206665,0.9503503327975676,0.0001454122746379124,0.9508179927088858,0.0001454480338635091,0.9512904450832298,0.0001454841778942922,0.9517679821256925,0.0001455207285088794,
0.9522508947756375,0.00014555770732492,0.9527394721549711,0.0001455951357527194,0.9532340009683443,0.0001456330349451717,0.9537347648542913,0.0001456714257439507,0.9542420436865237,0.0001457103286219321,
0.9547561128248467,0.0001457497636218381,0.9552772423154445,0.0001457897502911246,0.9558056960406009,0.0001458303076131549,0.9563417308182725,0.0001458714539347385,0.9568855954523308,0.0001459132068901443,
0.957437529734715,0.000145955583321735,0.9579977634012095,0.0001459985991974088,0.9585665150430652,0.0001460422695250747,0.9591439909772242,0.0001460866082644329,0.9597303840784822,0.0001461316282363772,
0.9603258725775256,3.429496445563442e-05,4.71390180589219e-05,0.02545140450924156,5.698005902228804e-05,0.09888482887297308,6.46665154872123e-05,0.1745672748051217,7.080319239456012e-05,0.2432554843126429,
7.580213079430914e-05,0.3035563500747286,7.994640198168628e-05,0.3560733630370389,8.343456150015368e-05,0.4018508725952348,8.640911033247278e-05,0.4419236498367694,8.897459090834537e-05,0.477195809384324,
9.120925845958284e-05,0.5084236299136395,9.317276496329196e-05,0.5362298756774806,9.491132926180562e-05,0.5611253276477577,9.646129005158931e-05,0.5835295355152852,9.785159643187924e-05,0.6037885335602649,
9.910558598615793e-05,0.6221892621956204,0.000100242275747997,0.6389710741614671,0.0001012773140638229,0.6543348452257979,0.0001022236924158808,0.6684501789683019,0.0001030922846881736,0.6814611164733789,
0.000103892260607767,0.6934906792512014,0.0001046314062227422,0.7046445019124451,0.0001053163748547551,0.7150137530267155,0.0001059528854651743,0.7246774971389701,0.0001065458808280527,0.7337046159216464,
0.0001070996546816144,0.7421553796866258,0.0001076179547185671,0.7500827400704289,0.000108104066598062,0.7575333991194183,0.0001085608829314993,0.7645486980661065,0.0001089909602824242,0.7711653599123298,
0.0001093965665385209,0.7774161128496166,0.0001097797204987463,0.7833302160485504,0.0001101422251266719,0.7889339050599856,0.0001104856956203599,0.7942507707081078,0.0001108115832166325,0.7993020827048136,
0.0001111211954666116,0.8041070671151008,0.0001114157135775683,0.8086831451312193,0.0001116962073042465,0.8130461392755103,0.0001119636477840458,0.8172104520762232,0.0001122189186395992,0.8211892213917328,
0.0001124628256154204,0.8249944558535597,0.0001126961049694284,0.8286371533240166,0.0001129194308029721,0.8321274047940588,0.0001131334214826889,0.8354744857604804,0.0001133386452827434,0.8386869368027751,
0.0001135356253556231,0.8417726348159535,0.0001137248441228562,0.8447388561361401,0.0001139067471630854,0.847592332612694,0.0001140817466633488,0.8503393015273689,0.0001142502244897405,0.8529855501323292,
0.0001144125349255261,0.8555364554704203,0.000114569007117969,0.8579970200494635,0.0001147199472693837,0.8603719038646763,0.000114865640603067,0.8626654531972963,0.0001150063531306337,0.8648817265611982,
0.0001151423332437714,0.8670245181211877,0.0001152738131504335,0.8690973788654197,0.0001154010101729216,0.8711036357789597,0.0001155241279231094,0.8730464092349987,0.0001156433573681646,0.8749286287938791,
0.000115758877798493,0.8767530475772949,0.0001158708577082165,0.8785222553652559,0.0001159794555972724,0.8802386905462134,0.0001160848207031594,0.8819046510357806,0.0001161870936694311,0.8835223042664072,
0.0001162864071572288,0.8850936963389427,0.0001163828864054445,0.8866207604170093,0.0001164766497444822,0.8881053244363083,0.0001165678090680486,0.8895491181932601,0.0001166564702669258,0.8909537798705544,
0.0001167427336282587,0.892320862051182,0.0001168266942035215,0.893651837267199,0.0001169084421479998,0.8949481031247646,0.0001169880630343355,0.896210987042823,0.0001170656381424237,0.8974417506390828,
0.0001171412447277278,0.8986415937936533,0.0001172149562698677,0.8998116584177496,0.0001172868427031662,0.9009530319522595,0.0001173569706306669,0.9020667506186153,0.0001174254035230018,0.9031538024423191,
0.0001174922019033529,0.9042151300675857,0.0001175574235196376,0.9052516333798847,0.0001176211235049476,0.9062641719516429,0.0001176833545271722,0.9072535673250101,0.0001177441669286597,0.9082206051443592,
0.0001178036088566899,0.9091660371500875,0.0001178617263854677,0.9100905830442862,0.0001179185636302828,0.910994932237941,0.0001179741628544263,0.9118797454885111,0.0001180285645694077,0.9127456564359945,
0.0001180818076289649,0.9135932730449129,0.0001181339293173255,0.9144231789590457,0.0001181849654321338,0.9152359347751857,0.0001182349503624295,0.9160320792416857,0.0001182839171620296,0.9168121303871106,
0.0001183318976186388,0.9175765865838867,0.0001183789223189884,0.9183259275514641,0.000118425020710278,0.9190606153031572,0.0001184702211581763,0.9197810950405108,0.000118514551001616,0.9204877959987503,
0.0001185580366045994,0.9211811322466059,0.0001186007034052168,0.9218615034435585,0.0001186425759620638,0.9225292955573307,0.0001186836779982295,0.9231848815442406,0.0001187240324430169,0.9238286219948466,
0.0001187636614715423,0.9244608657471405,0.0001188025865423529,0.9250819504693824,0.0001188408284331903,0.9256922032145273,0.0001188784072750181,0.9262919409480546,0.0001189153425844261,0.9268814710508878,
0.000118951653294512,0.9274610917989764,0.0001189873577843384,0.928031092821002,0.0001190224739070536,0.9285917555355758,0.0001190570190167592,0.9291433535691963,0.0001190910099942032,0.9296861531561591,
0.0001191244632713702,0.9302204135215227,0.0001191573948550362,0.9307463872481658,0.0001191898203493515,0.9312643206289019,0.0001192217549775109,0.9317744540045509,0.000119253213602564,0.9322770220888102,
0.0001192842107474197,0.9327722542807086,0.0001193147606140906,0.9332603749653779,0.0001193448771022213,0.9337416038038212,0.0001193745738269449,0.9342161560123183,0.0001194038641361023,0.9346842426320593,
0.0001194327611268635,0.9351460707895553,0.0001194612776617823,0.9356018439483411,0.0001194894263843151,0.936051762152442,0.0001195172197338335,0.9364960222620436,0.0001195446699601559,0.9369348181817692,
0.0001195717891376224,0.9373683410819361,0.0001195985891787351,0.9377967796131313,0.0001196250818473841,0.9382203201144163,0.0001196512787716769,0.9386391468154388,0.0001196771914563876,0.9390534420327067,
0.0001197028312950406,0.9394633863602434,0.0001197282095816411,0.9398691588548227,0.0001197533375220634,0.9402709372159495,0.0001197782262451063,0.9406688979607301,0.0001198028868132232,0.941063216593744,
0.0001198273302329317,0.9414540677720091,0.0001198515674649078,0.9418416254650982,0.000119875609433767,0.942226063110445,0.0001198994670375319,0.9426075537638465,0.0001199231511567857,0.9429862702451438,
0.0001199466726635097,0.943362385279036,0.0001199700424295984,0.9437360716309537,0.0001199932713350474,0.9441075022378914,0.0001200163702758047,0.9444768503340695,0.0001200393501712764,0.944844289571268,
0.0001200622219714742,0.9452099941336445,0.0001200849966637913,0.9455741388468206,0.0001201076852793909,0.9459368992809881,0.0001201302988991895,0.9462984518477588,0.0001201528486594163,0.9466589738904483,
0.000120175345756726,0.9470186437674541,0.0001201978014528419,0.9473776409283557,0.0001202202270787045,0.9477361459823309,0.0001202426340380965,0.9480943407584549,0.0001202650338107162,0.9484524083574091,
0.0001202874379546662,0.9488105331941014,0.0001203098581083251,0.9491689010306606,0.0001203323059915656,0.9495276989992437,0.0001203547934062821,0.9498871156140531,0.000120377332236188,0.9502473407719423,
0.0001203999344458418,0.950608565740949,0.0001204226120788594,0.9509709831360716,0.0001204453772552672,0.9513347868815795,0.0001204682421679518,0.951700172159119,0.0001204912190781571,0.9520673353408602,
0.0001205143203099825,0.9524364739069084,0.0001205375582438321,0.9528077863461863,0.000120560945308766,0.9531814720399889,0.0001205844939737027,0.9535577311273962,0.000120608216737424,0.9539367643517374,
0.0001206321261173331,0.9543187728872938,0.0001206562346369151,0.9547039581454467,0.0001206805548118562,0.9550925215594878,0.0001207050991347726,0.9554846643473366,0.0001207298800585086,0.9558805872514442,
0.0001207549099779616,0.9562804902552012,0.0001207802012103997,0.9566845722752253,0.0001208057659742362,0.9570930308289628,0.0001208316163662374,0.9575060616771147,0.0001208577643371385,0.9579238584404856,
0.0001208842216656548,0.9583466121909501,0.0001209109999308791,0.958774511016349,0.0001209381104830662,0.9592077395592462,0.0001209655644128141,0.9596464785296278,0.0001209933725186615,0.9600909041917725,
0.000121021545273132,0.9605411878256951,0.0001210500927872672,0.9609974951637525,0.0001210790247737047,0.9614599858031972,0.0001211083505083695,0.9619288125956789,0.0001211380787908609,0.9624041210149245,
0.0001211682179036361,0.9628860485040655,0.0001211987755701019,0.9633747238043349,0.0001212297589117476,0.9638702662671197};

//coeff 1 of final two stage procedure
Real ROCK2::fp1[46]={0.4102693550421609,0.3889624104727243,0.3804692420283886,0.3760815680865637,0.3735177579729938,0.3719340231904236,0.3708571145968057,0.3700947006022557,0.3695328931459086,0.3691085831661758,
0.368781324965233,0.3685244707068931,0.3683185599507446,0.3681542178682514,0.3680181997765286,0.3679084456991284,0.3678181571053212,0.3678314333608541,0.3677897070402892,0.3681800192470787,
0.3681272993461229,0.3680840569645587,0.3680522380648169,0.3680263578626069,0.3680061275157194,0.3679837719607466,0.3679668653311732,0.3679542340323301,0.367942933258425,0.3679349432021754,
0.3679290943359695,0.3679242023884676,0.3679207541681089,0.3679185472223537,0.367916869013064,0.3679158588043139,0.3679154592969145,0.3679154025286917,0.3679157536198652,0.3679163676763697,
0.3679171904021983,0.3679181786833088,0.3679192462983425,0.367920432307971,0.3679216942157868,0.3679229127010114};

//coeff 2 of final two stage procedure
Real ROCK2::fp2[46]={0.4495196112243335,0.4219428123056774,0.4084335547255627,0.4009301129475925,0.3963598727888637,0.3934034185789226,0.3913676516238603,0.3899091428928617,0.388827696299666,0.3880048656683555,
0.3873650613539532,0.3868583585730354,0.3864499054795832,0.3861178821587815,0.3858426294881124,0.3856144554520791,0.3854228843194507,0.3853156085078759,0.3851902798680153,0.3853705093720269,
0.3851957294861824,0.3850587241670235,0.3849515900397918,0.3848648995575697,0.38479450822313,0.38471172244074,0.384647872630907,0.3845979020112835,0.3845472492329918,0.3845088319509754,
0.3844789634065264,0.3844504092359686,0.3844285235163634,0.3844115660464609,0.3843957877945634,0.3843835767106759,0.3843727001444428,0.3843632394180673,0.3843553120908175,0.3843487746818896,
0.3843434800804404,0.3843392605995229,0.3843359163858929,0.3843331309176606,0.3843309056586355,0.3843292556220249};


Real DROCK2::recf2[184]={0.3005488243231985,0.03229561778843724,0.6512099103020813,0.6392740828523978,0.1881787048433069,0.1294908679614229,0.3692922899497312,0.7942384047165153,0.1266454719228368,0.2166410891457194,
0.2031748666628389,0.6897398524208641,0.09219088832162044,0.2919611374292846,0.1295377667266644,0.6475306974112423,0.07046743565977052,0.3554268168311739,0.090765655958888,0.6330845728770385,
0.05569721632404252,0.4082111635239982,0.0675873603555146,0.629870206918841,0.04523145105404241,0.4534919541952099,0.05267926402727696,0.6354404914695549,0.03750900894316795,0.4925487854409745,
0.0424187172939975,0.6449246344628944,0.03163592498422852,0.5265592512132928,0.0350084592084723,0.6562593108616718,0.02705675781155015,0.5563731440005397,0.02945185304996798,0.668208311194205,
0.02341331121971601,0.5826905173274041,0.02516196190460799,0.6801621264871356,0.02046467755479506,0.6060751338104287,0.02177157401321578,0.6918070079811171,0.01804385306038518,0.6269996819104969,
0.01904055637181734,0.7030168155266024,0.01603037918981687,0.6457867627448038,0.01680366913196333,0.7136405529152988,0.01433811555956923,0.662793143176811,0.01494746986542325,0.7237497904376623,
0.01290104693454032,0.6782069833369995,0.01338760195273448,0.7332517345924864,0.01167040637506997,0.6922580908837453,0.01206358679431103,0.7422067088089256,0.01059608600689598,0.7040977415114994,
0.01091371576681292,0.7491236879479691,0.009671969793341206,0.7157191213478938,0.009933289409588535,0.7568212787292203,0.008830880492643542,0.7253114860815413,0.009044704051776029,0.7625263511197745,
0.007497529753444489,0.7442970376537648,0.007648004349791438,0.7757801402936246,0.006446624996274999,0.7609681570976378,0.006555791376976517,0.7880046700358945,0.005602816140246905,0.775610944163273,
0.005683953042659415,0.7990879327137264,0.004914987115630155,0.7885753716314036,0.004976559648558256,0.8091578532120939,0.004346776360990677,0.8001306739546789,0.004394351755547853,0.8183261167651231,
0.003663188961270879,0.815282107874271,0.003696452253824382,0.8306098315072843,0.003129325162612305,0.8283070348886951,0.00315328263217451,0.8413983297008319,0.002704369633781579,0.8396225318494346,
0.002722063424297376,0.8509354578173889,0.002260423597593087,0.8525791014289577,0.002272630242688298,0.862034774036096,0.001917562864983125,0.8636040019677361,0.001926255663358357,0.8716262263454221,
0.001647248304620235,0.8730977440905955,0.001653606320243899,0.8799901295738884,0.001382826371800503,0.8832578837444272,0.001387265296404007,0.8890451857309281,0.001177354156527833,0.8919144694755784,
0.001180546848175413,0.8968431219167784,0.001014515810920781,0.8993770061880118,0.001016870670639793,0.9036250562189186,0.0008601033701745554,0.9070773083168101,0.0008617844954586736,0.9106798720695611,
0.0007384498387159457,0.913683985812873,0.0007396819312227429,0.9167779015832841,0.0006265364943395239,0.9202958616358844,0.0006274183910651552,0.9229216962591809,0.0005271988414776385,0.926712691547667,
0.0005278198525697253,0.9289229262885267,0.0004412803054223307,0.9327997241975136,0.0004417131675794807,0.934650367770951,0.0003683298903532936,0.9384778923722996,0.0003686300381467888,0.9400231060653412,
0.0003071781069080669,0.9437102190186861,0.0003073859624892617,0.9449993030326646,0.0002563479820628573,0.9484896107065751,0.0002564921728066633,0.9495657124903313,0.0002143136578947937,0.9528287843712028,
0.000214414083069791,0.9537286932946777,0.0001774960733422132,0.9570077589594689,0.000177564725268418,0.9577532838672067,0.0001461773410303843,0.9609306188294088,0.0001462237569103017,0.961544768057688,
0.0001212611744044664,0.9643727851516277,0.000121293027834231,0.9648823789107062};

Real DROCK2::recalph[46]={1.029217577908589,1.171142615787386,1.289178731889639,1.370039721728061,1.430279163156856,1.478209340802438,1.516334027924159,1.547598809712693,1.573716033722604,1.595931553195606,
1.6150939953884995,1.631813398460148,1.646519697230067,1.659614432610494,1.67129422100608,1.681847425139011,1.691412128447207,1.701553522119702,1.709840773771184,1.723292787981514,
1.736782707915586,1.748192229130831,1.758055997758203,1.766639784928805,1.774194682070617,1.783976225730842,1.792257041217689,1.799364729395389,1.807425652081097,1.814205683210333,
1.819986883764463,1.826122183655859,1.83130311714655,1.835740269183132,1.840289126639967,1.844167818112841,1.848031285616063,1.85175987575436,1.8552799194745,1.858548029126852,
1.861546773456168,1.864275803993039,1.866744876869206,1.869115278174705,1.871333924176762,1.873275509706936};



