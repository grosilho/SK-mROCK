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