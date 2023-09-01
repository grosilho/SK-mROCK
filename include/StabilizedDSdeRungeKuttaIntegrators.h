#ifndef STABILIZEDDSDERUNGEKUTTAINTEGRATORS_H
#define STABILIZEDDSDERUNGEKUTTAINTEGRATORS_H

#include "DSdeRungeKuttaIntegrator.h"
#include "StabilizedOdeRungeKuttaIntegrators.h"

class SKROCK: public DSdeRungeKuttaIntegrator, 
              public RKC1
{
public:
    SKROCK(Parameters* param_, DSde* sde_);
    virtual ~SKROCK();
    
    void step(const Real t, const Real& h);
    void update_n_stages_and_h(Real& h);
};

class SKmROCK: public virtual MultirateDSdeRungeKuttaIntegrator, public virtual mRKC
{
public:
    SKmROCK(Parameters* param_, MultirateDSde* msde_);
    virtual ~SKmROCK();

protected:    
    void g_eta(Real t, Vector& x, Vector& gx);
    void step(const Real t, const Real& h);    
    void update_n_stages_and_h(Real& h);      
    virtual void disp_step_info(Real& t, Real& h);

};

#endif /* STABILIZEDDSDERUNGEKUTTAINTEGRATORS_H */

