#ifndef FASTSLOWODERUNGEKUTTAINTEGRATOR_H
#define FASTSLOWODERUNGEKUTTAINTEGRATOR_H

#include "MainHeader.h"
#include "OdeRungeKuttaIntegrator.h"

class MultirateOdeRungeKuttaIntegrator: public OdeRungeKuttaIntegrator
{
public:
    MultirateOdeRungeKuttaIntegrator(Parameters* param_, MultirateOde* mode_);
    virtual ~MultirateOdeRungeKuttaIntegrator();
                      
    virtual void print_integration_info();   
        
protected:
    virtual bool update_rho(Real t);      //updates spectral radius
    virtual unsigned int rho(Real t, int& iter_F, int& iter_S);//power method approximating rho        
  
    virtual void reinit_statistics();
    
    virtual void disp_step_info(Real& t, Real& h, bool accepted);
    
protected:
    MultirateOde* mode;
    
    //equation variables
    Vector* eigenvector_F;
    Vector* eigenvector_S;
//    Vector* integr_add[3]; //working vectors
    
    //integration statistics
    int max_rho_F;             ///<Maximal spectral radius computed.
    int min_rho_F;             ///<Minimal spectral radius computed.
    int max_rho_S;             ///<Maximal spectral radius computed.
    int min_rho_S;             ///<Minimal spectral radius computed.
    int m_max;                ///<Maximal number of stages used.   
    int m_avg;                ///<Average number of stages used.   
    int n_fF_eval_rho;        ///<Number of righ hand side evaluations for the spectral radius computation in the power method.
    int n_fF_eval;            ///<Number of righ hand side evaluations for time integration.
    int n_fS_eval_rho;        ///<Number of righ hand side evaluations for the spectral radius computation in the power method.
    int n_fS_eval;            ///<Number of righ hand side evaluations for time integration.
   
    //step variables
    Real eigmax_F;       ///<Spectral radius of current time step.    
    Real eigmax_S;       ///<Spectral radius of current time step.    
    int m;               
};

#endif /* FASTSLOWODERUNGEKUTTAINTEGRATOR_H */

