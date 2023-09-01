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

#endif /* STABILIZEDDSDERUNGEKUTTAINTEGRATORS_H */

