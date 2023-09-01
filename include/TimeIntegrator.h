#ifndef TIMEINTEGRATOR_H
#define	TIMEINTEGRATOR_H

#include "Parameters.h"
#include "Ode.h"


class TimeIntegrator
{
public:
    TimeIntegrator(Parameters* param_, Ode* ode_);
    virtual ~TimeIntegrator();
 
    virtual bool integrate() =0; //compute next yn 
    virtual void convergence_test() =0;
    virtual Vector solution() =0;
    
    virtual void print_integration_info() =0;

protected:    
    Real get_cpu_time();
    void output_solution(Real t, Vector* y);
    void output_final_solution(Vector* y);
    void read_reference_solution(Vector* refsol);
    
    virtual void compute_errors();

    Parameters* param;
    Ode* ode;    
    Real elapsed_time;
    
    unsigned int n_output;
};

#endif	/* TIMEINTEGRATOR_H */

