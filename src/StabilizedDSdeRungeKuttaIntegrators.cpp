#include "StabilizedDSdeRungeKuttaIntegrators.h"
#include "ChebyshevMethods.h"

SKROCK::SKROCK(Parameters* param_, DSde* sde_)
:OdeRungeKuttaIntegrator(param_,sde_),
 DSdeRungeKuttaIntegrator(param_,sde_),
 RKC1(param_,sde_)
{
    reinit_integrator();
}

SKROCK::~SKROCK()
{
}

void SKROCK::update_n_stages_and_h(Real& h)
{
    RKC1::update_n_stages_and_h(h);
}


void SKROCK::step(const Real t, const Real& h)
{
    Vector*& Q = integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> nu;
    vector<Real> mu;
    vector<Real> kappa;
    
    ChebyshevMethods::CoefficientsSKROCK(mu, nu, kappa, c, s, damping);
    
    if(sde->get_noise_type()==DIAGONAL)
    {
        sde->g(t,*yn,*Kjm2);   
        *Q = Kjm2->cwiseProduct(I->get_Ir());
    }
    else
    {
        sde->g(t,*yn,G);
        *Q = G*I->get_Ir();
    }
    
    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    *Kjm2 = *yn + nu[0]*(*Q);
    sde->f(t,*Kjm2,*Kj); 
    *Kj *= mu[0]*h;
    *Kj += kappa[0]*(*Q);
    *Kj += *Kjm1;
    
    // Stages j=2,...,s
    for(unsigned int j=2;j<=s;j++)
    {
        swap_ptr=Kjm2;
        Kjm2=Kjm1;
        Kjm1=Kj;
        Kj=swap_ptr;               
                
        sde->f(t+h*c[j-2],*Kjm1,*Kj);           
        *Kj *= h*mu[j-1];
        *Kj += nu[j-1]*(*Kjm1);
        *Kj += kappa[j-1]*(*Kjm2);                                   
    }   
    
    swap_ptr=ynpu;
    ynpu=Kj;
    Kj=swap_ptr;    
    
    n_f_eval += s;
    n_g_eval ++;
}


SKmROCK::SKmROCK(Parameters* param_, MultirateDSde* msde_)
:OdeRungeKuttaIntegrator(param_,msde_),
 MultirateOdeRungeKuttaIntegrator(param_,msde_),
 DSdeRungeKuttaIntegrator(param_, msde_), 
 MultirateDSdeRungeKuttaIntegrator(param_,msde_),
 mRKC(param_,msde_)
{
    reinit_integrator();
}

SKmROCK::~SKmROCK()
{
}

void SKmROCK::step(const Real t, const Real& h)
{
    Vector*& Q = integr[0];
    Vector*& Kjm1= integr[1];
    Vector*& Kjm2= integr[2];
    Vector*& Kj= integr[3];
    Vector* swap_ptr=0;   
    
    vector<Real> c;
    vector<Real> nu;
    vector<Real> mu;
    vector<Real> kappa;
    
    ChebyshevMethods::CoefficientsSKROCK(mu, nu, kappa, c, s, damping);

    g_eta(t,*yn,*Q);

    // Stage 0
    *Kjm1 = *yn;
        
    // Stage 1
    *Kjm2 = *yn + nu[0]*(*Q);
    // f_eta(t,*Kjm2,*fn); 
    // *Kj = mu[0]*h*(*fn);
    // *Kj += kappa[0]*(*Q);
    // *Kj += *Kjm1;

    f_eta(t,*Kjm2,*Kj); 
    *Kj *= mu[0]*h;
    *Kj += kappa[0]*(*Q);
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
    
    n_fS_eval += s;
    n_fF_eval += s*m;
}

void SKmROCK::g_eta(Real t, Vector& x, Vector& Qeta)
{
    unsigned r = m/2;
    Vector*& vjm2 = integr[4];
    Vector*& vjm1 = integr[5];
    Vector*& vj = integr[6];
    Vector*& gdW = integr[7];
    Vector* swap_ptr=0;   

    vector<Real> c;
    vector<Real> beta;
    vector<Real> alpha;
    vector<Real> gamma;
    Real theta;

    ChebyshevMethods::CoefficientsRKC1mod(alpha, beta, gamma, theta, c, m, damping);

    if(sde->get_noise_type()==DIAGONAL)
    {
        sde->g(t,x,*vjm2);   
        *gdW = vjm2->cwiseProduct(I->get_Ir());
    }
    else
    {
        sde->g(t,x,G);
        *gdW = G*I->get_Ir();
    }

    *vjm1 = x;
    *vjm2 = x + beta[0]*theta*eta*(*gdW);
    msde->fF(t,*vjm2,*vj);
    *vj *= alpha[0]*eta;
    *vj += gamma[0]*theta*eta*(*gdW);
    *vj += x;
    for(unsigned j=2;j<=r;j++)
    {
        swap_ptr=vjm2;
        vjm2=vjm1;
        vjm1=vj;
        vj=swap_ptr;   

        msde->fF(t,*vjm1,*vj);
        *vj *= alpha[j-1]*eta;
        *vj += beta[j-1]*(*vjm1);
        *vj += gamma[j-1]*(*vjm2);
    }
    Qeta = *vj;

    *vjm1 = x;
    msde->fF(t,*vjm1,*vj);
    *vj *= alpha[0]*eta;
    *vj += x;
    for(unsigned j=2;j<=r;j++)
    {
        swap_ptr=vjm2;
        vjm2=vjm1;
        vjm1=vj;
        vj=swap_ptr;   

        msde->fF(t,*vjm1,*vj);
        *vj *= alpha[j-1]*eta;
        *vj += beta[j-1]*(*vjm1);
        *vj += gamma[j-1]*(*vjm2);
    }

    Qeta -= *vj;
    Qeta *= (1./eta);
}
    
void SKmROCK::update_n_stages_and_h(Real& h)
{
    /**
     * Updates the number of stages needed by SKmROCK. 
     */
    
    unsigned int stages_limit = 1e6;
    
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
    }

    m = ceil(sqrt(6.*eigmax_F*h/beta/beta/s/s+1.))+safe_add;
    m = m + m%2;
    eta = 6.*m*m*h/(m*m-1.)/beta/s/s;
    
//    write_parameters(h);
    
    s_max= max(s_max,s);
    m_max= max(m_max,m);
    s_avg += s;
    m_avg += m;
}
    
void SKmROCK::disp_step_info(Real& t, Real& h)
{
    MultirateDSdeRungeKuttaIntegrator::disp_step_info(t,h);
}