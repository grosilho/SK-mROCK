#ifndef ODERUNGEKUTTAINTEGRATOR_H
#define	ODERUNGEKUTTAINTEGRATOR_H

#include "MainHeader.h"
#include "TimeIntegrator.h"

class OdeRungeKuttaIntegrator: public TimeIntegrator
{
public:
    OdeRungeKuttaIntegrator(Parameters* param_, Ode* ode_);
    virtual ~OdeRungeKuttaIntegrator();
                  
    virtual bool integrate(); //solves problem
    virtual bool integrate(const Vector& y0, Real t0, Real T, Real dt, bool new_integration); //solves problem with given data
    virtual void convergence_test();
    Vector solution();
    
    virtual void reinit_integrator();
    
    virtual void print_integration_info();
    
    unsigned int get_n_f_eval();
    unsigned int get_n_f_eval_rho_or_Nit();
    unsigned int get_n_steps();
      
protected:
    virtual void step(const Real t, const Real& h) = 0;
    
    virtual bool update_rho(Real t);      //updates spectral radius
    virtual void update_n_stages_and_h(Real& h) = 0; //given rho, computes number of stages
    virtual unsigned int  rho(Real t, int& iter);//power method approximating rho
    virtual void reinit_statistics();
    
    virtual void disp_step_info(Real& t, Real& h);
    
protected:
    //equation variables
    Vector* yn;
    Vector* ynpu;
    Vector* eigenvector;
    Vector* integr[5]; //working vectors
    Real t;
    Real tend;
    Real h;
        
    //integration settings
    bool internal_rho;   ///<If true the ODE's spectral radius is computed by the internal nonlinear power method. Else a method rho is provived by the ODE class.
    unsigned int rho_freq;  ///<Computes the spectral radius at least onece every rho_freq steps.
    bool verbose;        ///<If true prints info about the timestep.
    int output_frequency;
    
    //integration statistics
    int max_rho;             ///<Maximal spectral radius computed.
    int min_rho;             ///<Minimal spectral radius computed.
    int s_max;               ///<Maximal number of stages used.
    Real s_avg;              ///<Average number of stages used, per step.
    int n_f_eval_rho;        ///<Number of righ hand side evaluations for the spectral radius computation in the power method or in the Newton algorithm of implicit schemes.
    int n_f_eval;            ///<Number of righ hand side evaluations for time integration.
    int n_steps;             ///<Number of time steps.
    
    //step variables
    Real eigmax;         ///<Spectral radius of current time step.
    int s;               ///<Number of stages for RKC. Number of stages minus 2 for ROCK2.
    int sold;            ///<Previous number of stages.
    int nrho;            ///<Number of accpeted steps after last spectral radius estimation.
    bool last;           ///<Is true if the current step is the last one.
    bool Astable;        ///<Tells if the RK method is A stable.
    
    bool rho_outdated;
};



#endif	/* ODERUNGEKUTTAINTEGRATOR_H */