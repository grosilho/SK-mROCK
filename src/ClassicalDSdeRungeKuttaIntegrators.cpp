#include "ClassicalDSdeRungeKuttaIntegrators.h"

EulerMaruyama::EulerMaruyama(Parameters* param_, DSde* sde_)
:OdeRungeKuttaIntegrator(param_,sde_),
 DSdeRungeKuttaIntegrator(param_,sde_),
 ExplicitEuler(param_,sde_)
{
}

EulerMaruyama::~EulerMaruyama()
{
}

void EulerMaruyama::update_n_stages_and_h(Real& h)
{
    // This is WRONG. As EM has a stronger stability condition than EE
    // But it would need the computation of the spectral radius of the diffusion term
    ExplicitEuler::update_n_stages_and_h(h);
}

void EulerMaruyama::step(const Real t, const Real& h)
{
    Vector*& k1= integr[0];
        
    sde->f(t,*yn,*k1);
    *ynpu = *yn + h*(*k1);

    if(sde->get_noise_type()==DIAGONAL)
    {
        sde->g(t,*yn,*k1);   
        *ynpu += k1->cwiseProduct(I->get_Ir());
    }
    else
    {
        sde->g(t,*yn,G);
        *ynpu += G*I->get_Ir();
    }
    
    n_f_eval ++;
    n_g_eval ++;
}

PlatenScheme::PlatenScheme(Parameters* param_, DSde* sde_)
:OdeRungeKuttaIntegrator(param_,sde_),
 DSdeRungeKuttaIntegrator(param_,sde_),
 ExplicitEuler(param_,sde_)
{
    needDoubleIntegral = true;
}

PlatenScheme::~PlatenScheme()
{
}

void PlatenScheme::update_n_stages_and_h(Real& h)
{
    // This is WRONG. But I dont want to estimate the specral radius of the diffusion term
    ExplicitEuler::update_n_stages_and_h(h);
}

void PlatenScheme::step(const Real t, const Real& h)
{
    Vector*& k1= integr[0];
//    Vector*& k2= integr[1];
    Vector*& kk= integr[2];
    Vector*& gj= integr[3];
    Vector*& GI= integr[4];
    
    Real sqrh = sqrt(h);
    
    sde->f(t,*yn,*k1);
    *k1 *= h;
    *k1 += *yn;

    if(sde->get_noise_type()==DIAGONAL)
    {
        sde->g(t,*yn,*GI);
        *ynpu = *k1 + GI->cwiseProduct(I->get_Ir());
        
        *k1 = *yn + sqrh*(*GI);
        sde->g(t,*k1,*gj);
        *gj -= *GI;
        
        *k1 = I->get_Ir().array().pow(2);
        *k1 = (k1->array()-h)/2./sqrh;
        
        *ynpu += gj->cwiseProduct(*k1);
        
        n_g_eval++;
    }
    else if(sde->get_noise_type()==COMMUTATIVE)
    {
        sde->g(t,*yn,G);
        *GI = G*I->get_Ir();
        *ynpu = *k1 + *GI;

        sde->g(t,*ynpu,G);
        *ynpu += (G*I->get_Ir()-*GI)/2.;
        
        n_g_eval += 2*sde->brownian_size();
    }
    else// GENERAL noise
    {
        sde->g(t,*yn,G);
        *ynpu = *k1 + G*I->get_Ir();
        
//        cout<<"norm I "<<I->get_Ir().norm()<<endl;
//        cout<<"norm ynpu "<<ynpu->norm()<<endl;

        const unsigned int m = sde->brownian_size();
        const unsigned int s_size = sde->get_system_size();
        

        for(unsigned int k=1;k<=m;k++)
        {
            *kk = *k1 + G.block(0,k-1,s_size,1)*sqrh;
            
//            cout<<"norm kk "<<kk->norm()<<endl;
            
            for(unsigned int j=1;j<=m;j++)
            {
                sde->g(t,*kk,*gj,j);
                *ynpu += (*gj-G.block(0,j-1,s_size,1))*I->get_Ipq()(k,j)/sqrh;
//                cout<<"norm Ipq "<<I->get_Ipq()(k,j)<<endl;
            }
        }
        n_g_eval += m*m;
    }
    
    n_f_eval ++;
}