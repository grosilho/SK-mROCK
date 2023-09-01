#ifndef CLASSICALDSDERUNGEKUTTAINTEGRATORS_H
#define CLASSICALDSDERUNGEKUTTAINTEGRATORS_H

#include "DSdeRungeKuttaIntegrator.h"
#include "ClassicalOdeRungeKuttaIntegrators.h"

class EulerMaruyama: public DSdeRungeKuttaIntegrator, 
                     public ExplicitEuler
{
public:
    EulerMaruyama(Parameters* param_, DSde* sde_);
    virtual ~EulerMaruyama();
    
    void step(const Real t, const Real& h);
    void update_n_stages_and_h(Real& h);
};


class PlatenScheme: public DSdeRungeKuttaIntegrator,
                    public ExplicitEuler
{
public:
    PlatenScheme(Parameters* param_, DSde* sde_);
    virtual ~PlatenScheme();
    
    void step(const Real t, const Real& h);
    void update_n_stages_and_h(Real& h);
};


#endif /* CLASSICALDSDERUNGEKUTTAINTEGRATORS_H */

