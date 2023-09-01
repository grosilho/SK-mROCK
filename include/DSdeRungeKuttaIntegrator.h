#ifndef DSDERUNGEKUTTAINTEGRATOR_H
#define DSDERUNGEKUTTAINTEGRATOR_H

#include "OdeRungeKuttaIntegrator.h"
#include "DSde.h"
#include "StochasticIntegrals.h"


class DSdeRungeKuttaIntegrator: public virtual OdeRungeKuttaIntegrator
{
public:
    DSdeRungeKuttaIntegrator(Parameters* param_, DSde* sde_);
    virtual ~DSdeRungeKuttaIntegrator();
                  
    virtual bool integrate(); //solves problem along a random path
    virtual bool integrate(HistoryStochasticIntegrals& hsi); //solves problem along a given path
    virtual void convergence_test();
    
    virtual void print_integration_info();
    
    bool need_double_integral();
    
protected:
    virtual void step(const Real t, const Real& h) = 0;
    virtual void update_n_stages_and_h(Real& h) = 0;
    
    virtual void reinit_statistics();
    virtual void disp_step_info(Real& t, Real& h);
    
    DSde* sde;
    StochasticIntegral* I;
    Matrix G;
    int n_g_eval;
    bool needDoubleIntegral;
};

#endif /* DSDERUNGEKUTTAINTEGRATOR_H */

