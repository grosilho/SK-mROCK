#include <fstream>

#include "MultirateOdeRungeKuttaIntegrator.h"

MultirateOdeRungeKuttaIntegrator::MultirateOdeRungeKuttaIntegrator(Parameters* param_, MultirateOde* mode_)
:OdeRungeKuttaIntegrator(param_, mode_), mode(mode_)
{
    eigenvector_F = new Vector(mode->get_system_size());
    eigenvector_S = new Vector(mode->get_system_size());
    
//    reinit_integrator();
}

MultirateOdeRungeKuttaIntegrator::~MultirateOdeRungeKuttaIntegrator()
{
    delete eigenvector_F;
    delete eigenvector_S;
}

void MultirateOdeRungeKuttaIntegrator::print_integration_info()
{
    string rho =u8"\u03C1";
    string delta = u8"\u0394";
    
    s_avg /= n_steps;
    m_avg /= n_steps;
    
    cout<<"\n-------------------   Integration Info   -------------------"<<endl;
    
    if(Astable)
        cout<<"The spectral radius has not been computed, "<<param->rk_name<<" is A-stable."<<endl;
    else
    {
        cout<<"Max "<<rho<<"F: "<<max_rho_F<<endl;
        cout<<"Min "<<rho<<"F: "<<min_rho_F<<endl;
        cout<<"Max "<<rho<<"S: "<<max_rho_S<<endl;
        cout<<"Min "<<rho<<"S: "<<min_rho_S<<endl;
        cout<<"Number of fF eval. for "<<rho<<"F: "<<n_fF_eval_rho<<endl;
        cout<<"Number of fS eval. for "<<rho<<"S: "<<n_fS_eval_rho<<endl;
    }
    cout<<"Max s: "<<s_max<<endl;
    cout<<"Mean s: "<<s_avg<<endl;
    cout<<"Max m: "<<m_max<<endl;
    cout<<"Mean m: "<<m_avg<<endl;
    cout<<"fF evaluations = "<<n_fF_eval<<endl;
    cout<<"fS evaluations = "<<n_fS_eval<<endl;
    cout<<"Maximal "<<delta<<"t used: "<<dt_max<<endl;
    cout<<"Number of steps: "<<n_steps<<endl;
    cout<<"Accepted steps: "<<acc_steps<<endl;
    cout<<"Rejected steps: "<<rej_steps<<endl; 
    cout<<"Elapsed time: "<<elapsed_time<<endl;
    cout<<"------------------------------------------------------------\n"<<endl;
}

void MultirateOdeRungeKuttaIntegrator::reinit_statistics()
{    
    OdeRungeKuttaIntegrator::reinit_statistics();
    
    m=0;
    m_max=0;
    m_avg=0.;
    
    max_rho_F=numeric_limits<int>::min();
    min_rho_F=numeric_limits<int>::max();
    max_rho_S=numeric_limits<int>::min();
    min_rho_S=numeric_limits<int>::max();
            
    n_fF_eval_rho=0;
    n_fF_eval=0;
    n_fS_eval_rho=0;
    n_fS_eval=0;
}

bool MultirateOdeRungeKuttaIntegrator::update_rho(Real t)
{
    /**
     * A new spectral radius is computed. Either with the estimation given by
     * ODE::rho or with the RungeKuttaIntegrator::rho internal power method.
     */
    
    if(Astable || (mode->is_rho_constant() && n_steps>0))
        return true;

    if(verbose) cout<<"\n--------------   Spectral Radius Estimation   --------------"<<endl;

    int iter_F=0, iter_S=0;
    unsigned int rho_conv;
    //Computed externally by MultirateODE::rho
    if (!internal_rho)
        mode->rho(t,*yn,eigmax_F,eigmax_S);
    //Computed internally by this->rho
    else
    {
        rho_conv = this->rho(t, iter_F, iter_S);

        if(rho_conv<=1)
        {
            cout<<"ERROR: convergence failure in spectral radius computation. "
                    "Fail computing "<<(rho_conv==0 ? "rho_F":"rho_S")<<endl;
            return false;
        }
        else if(rho_conv==2)
            cout<<"WARNING: augment number of iterations in spectral radius computation."<<endl;
    }

    eigmax = eigmax_S+eigmax_F; //usually smaller

    max_rho_F = max(max_rho_F,(int)eigmax_F+1);
    max_rho_S = max(max_rho_S,(int)eigmax_S+1);
    min_rho_F = min(min_rho_F,(int)eigmax_F+1);
    min_rho_S = min(min_rho_S,(int)eigmax_S+1);

    nrho=0;

    if(verbose) 
        cout<<scientific<<"Spectral radius estimation (F, S): ("<<eigmax_F<<", "<<eigmax_S<<")"<<endl;
    if(internal_rho && verbose)
        cout<<"Power method converged in ("<<iter_F<<", "<<iter_S<<") iterations."<<endl;

    
    if(verbose) cout<<"------------------------------------------------------------\n"<<endl;
    
    rho_outdated=false;
    
    return true;
}

unsigned int MultirateOdeRungeKuttaIntegrator::rho(Real t, int& iter_F, int& iter_S)
{
    /**
    *     rho computes eigmax, a close upper bound of the
    *     spectral radius of the Jacobian matrix using a 
    *     power method (J.N. Franklin (matrix theory)). 
    *     The algorithm used is a small change (initial vector
    *     and stopping criteria) of that of
    *     Sommeijer-Shampine-Verwer, implemented in RKC.
    */

    Real eigmaxo,sqrtu,znor,ynor,quot,dzyn,dfzfn;
    
    const int maxiter=100;
    const Real safe=1.05;
    const Real tol = 1e-3;
    
    sqrtu= sqrt(uround);

// ------ The initial vectors for the power method are yn --------
//       and yn+c*f(v_n), where vn=f(yn) a perturbation of yn 
//       (if n_steps=0) or a perturbation of the last computed
//       eigenvector (if n_steps!=0). 
    
    Vector*& fn = integr[0];
    Vector* z  = integr[1];
    
    srand(time(NULL));
    
    for(unsigned int fs=0;fs<=1;fs++)
    {
        Vector* eigenvectorFS = (fs==0 ? eigenvector_F:eigenvector_S);        
        
        if(n_steps==0)//take a random vector as first vector of power method
        {
            for(int i=0;i<mode->get_system_size();i++)
                (*z)(i) = ((rand()%2)-0.5)*((double)rand())/((double)RAND_MAX);
        }
        else//take the last vector of last power method call
            *z = *eigenvectorFS;

        if(fs==0)
            mode->fF(t,*yn,*fn);
        else
            mode->fS(t,*yn,*fn);

        // ------ Perturbation.--------
        ynor= yn->norm();
        znor= z->norm();

        // Building the vector z so that the difference z-yn is small
        if(ynor!=0.0 && znor!=0.0)
        {
            dzyn=ynor*sqrtu;
            quot=dzyn/znor;
            (*z) *= quot;
            (*z) += *yn;
        }
        else if(ynor!=0.0)
        {
            dzyn=ynor*sqrtu;
            *z=*yn;
            (*z) *= 1.+sqrtu;
        }
        else if(znor!=0.0)
        {
            dzyn=sqrtu;
            quot=dzyn/znor;
            (*z) *= quot;
        }
        else
        {
            dzyn=sqrtu*sqrt(z->size());
            for(int i=0;i<mode->get_system_size();i++)
                (*z)(i) += sqrtu;
        }
        //here dzyn=||z-yn|| and z=yn+(small perturbation)
        //dzyn=||z-yn|| will be always true, even with new z in the loop
        //Start the power method for non linear operator rhs

    //    eigmax=0.0;
        
        int& iter = (fs==0 ? iter_F:iter_S);
        int& n_g_eval_rho = (fs==0 ? n_fF_eval_rho:n_fS_eval_rho);
        Real& eigmaxFS = (fs==0 ? eigmax_F:eigmax_S);
        
        for(iter=1;iter<=maxiter;iter++)
        {
            if(fs==0)
                mode->fF(t,*z,*eigenvectorFS);
            else
                mode->fS(t,*z,*eigenvectorFS);
            n_g_eval_rho++;

            (*eigenvectorFS) -= *fn; //dz is the new perturbation, not normalized yet
            dfzfn= eigenvectorFS->norm();
    //        dfzfn= norm(*eigenvectorFS);
            
            eigmaxo=eigmaxFS;
            eigmaxFS=dfzfn/dzyn; //approximation of the Rayleigh quotient (not with dot product but just norms)
            eigmaxFS=safe*eigmaxFS;            

            if(abs(eigmaxFS-eigmaxo)<= eigmaxFS*tol)
            {
                //The last perturbation is stored. It will very likely be a
                // good starting point for the next rho call.
                *eigenvectorFS=*z;
                (*eigenvectorFS) -= *yn;    
                break;
            }
            if (dfzfn==0.0)
            {
                cout<<"ERROR: Convergence failure in the spectral radius computation."<<endl;
                return fs;
            }

            quot=dzyn/dfzfn;
            *z=*eigenvectorFS;
            (*z) *= quot;
            (*z) += *yn; //z is built so that dzyn=||z-yn|| still true

        }
    }
    
    if(iter_F==maxiter+1 || iter_S==maxiter+1)
        return 2;
    else
        return 3;
    
}

void MultirateOdeRungeKuttaIntegrator::disp_step_info(Real& t, Real& h, bool accepted)
{
    cout<<setprecision(4)<<scientific;
    
    string delta = u8"\u0394";
    
    cout << scientific;
    
    cout<<"Step t = "<<setw(6)<<setprecision(4)<<t<<", "<<delta<<"t = "<<setw(8)<<setprecision(6)<<h
    <<", s = "<<setw(3)<<s<<", m = "<<setw(3)<<m
    <<" and |y_n+1| = "<<setw(7)<<setprecision(4)
    <<ynpu->lpNorm<Eigen::Infinity>()<<". ";
    
    if(accepted)
        cout<<" Accepted "; 
    else
        cout<<" Rejected ";
    err_control.disp_info();
    cout<<endl;
}