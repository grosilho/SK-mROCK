#ifndef CLASSICALODERUNGEKUTTAINTEGRATORS_H
#define CLASSICALODERUNGEKUTTAINTEGRATORS_H

#include "OdeRungeKuttaIntegrator.h"

class ExplicitEuler : public virtual OdeRungeKuttaIntegrator
{
public:
    ExplicitEuler(Parameters *param_, Ode *ode_);
    virtual ~ExplicitEuler();

    virtual void step(const Real t, const Real &h);
    virtual void update_n_stages_and_h(Real &h);
};

#endif /* CLASSICALODERUNGEKUTTAINTEGRATORS_H */
