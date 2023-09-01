#ifndef ODEPROBLEMS_H
#define	ODEPROBLEMS_H

#include "Ode.h"

class DahlquistTestProblem: public virtual Ode
{
public:
    DahlquistTestProblem();
    virtual ~DahlquistTestProblem(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    void AN_df(Real t, Vector& x, Matrix& dfx); 
        
protected:
     static Real lambda;
     static Real xi;
};

class ScalarNonStiffNonLinearTest: public virtual Ode
{
public:
    ScalarNonStiffNonLinearTest();
    virtual ~ScalarNonStiffNonLinearTest(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    void AN_df(Real t, Vector& x, Matrix& dfx); 
};

class NeuronCable: public virtual Ode
{
public:
    NeuronCable();
    virtual ~NeuronCable(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    
protected:
    static Real nu;
    static Real beta;
};

class Brusselator: public virtual Ode
{
public:
    Brusselator();
    virtual ~Brusselator(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
        
protected:
    static Real alpha;
};

class PDEBrusselator: public virtual Ode
{
public:
    PDEBrusselator();
    virtual ~PDEBrusselator(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
            
protected:
    static Real A;
    static Real B;
    static Real alpha;
    static unsigned int Nu;
    static unsigned int Nv;
    static Real Hu;
    static Real Hv;
};

class PopulationDynamics: public virtual Ode
{
public:
    PopulationDynamics();
    virtual ~PopulationDynamics(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
            
protected:
    Real lambda1;
    Real lambda2;
    Real alpha;
};


#endif	/* ODEPROBLEMS_H */