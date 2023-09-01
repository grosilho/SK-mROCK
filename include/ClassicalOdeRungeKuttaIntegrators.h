#ifndef CLASSICALODERUNGEKUTTAINTEGRATORS_H
#define CLASSICALODERUNGEKUTTAINTEGRATORS_H

#include "OdeRungeKuttaIntegrator.h"


class ExplicitEuler: public virtual OdeRungeKuttaIntegrator
{
public:
    ExplicitEuler(Parameters* param_, Ode* ode_);
    virtual ~ExplicitEuler();
    
    virtual void step(const Real t, const Real& h);
    virtual void update_n_stages_and_h(Real& h);
};

class ImplicitEuler: public OdeRungeKuttaIntegrator
{
public:
    ImplicitEuler(Parameters* param_, Ode* ode_);
    virtual ~ImplicitEuler();
    
    virtual void step(const Real t, const Real& h);
    void update_n_stages_and_h(Real& h);
    
    virtual void disp_step_info(Real& t, Real& h, bool accepted);
    
private:
    Real Newton_tol;
    unsigned int Newton_max_iter;
    unsigned int Newton_iter;
    unsigned int lin_solv_iter;
    Matrix J;
    Matrix I;
    SpMatrix spJ;
    SpMatrix spI;
};

class ExplicitMidpoint: public OdeRungeKuttaIntegrator
{
public:
    ExplicitMidpoint(Parameters* param_, Ode* ode_);
    virtual ~ExplicitMidpoint();
    
    virtual void step(const Real t, const Real& h);
    void update_n_stages_and_h(Real& h);
};

class ImplicitMidpoint: public OdeRungeKuttaIntegrator
{
public:
    ImplicitMidpoint(Parameters* param_, Ode* ode_);
    virtual ~ImplicitMidpoint();
    
    virtual void step(const Real t, const Real& h);
    void update_n_stages_and_h(Real& h);
  
private:
    Real Newton_tol;
    unsigned int Newton_max_iter;
    unsigned int lin_solv_iter;
    Matrix J;
    Matrix I;
    SpMatrix spJ;
    SpMatrix spI;
};

class RungeKutta4: public OdeRungeKuttaIntegrator
{
public:
    RungeKutta4(Parameters* param_, Ode* ode_);
    virtual ~RungeKutta4();
    
    virtual void step(const Real t, const Real& h);
    void update_n_stages_and_h(Real& h);
};

#endif /* CLASSICALODERUNGEKUTTAINTEGRATORS_H */

