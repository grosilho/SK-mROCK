#ifndef PARAREAL_H
#define PARAREAL_H

#include "Parameters.h"
#include "Ode.h"
#include "OdeRungeKuttaIntegrator.h"

class Parareal
{
public:
    Parareal(Parameters param_);
    virtual ~Parareal();
    
    void integrate();
    
protected:
    void write_outer_y(Ode* ode, unsigned int iter);
    
    Real get_cpu_time();
    
    Parameters param;
    
    vector<Vector> outer_y;    
    vector<Vector> inner_y;
    vector<Vector> outer_int_sol;
    vector<Real> outer_t;
    
    
};

#endif /* PARAREAL_H */

