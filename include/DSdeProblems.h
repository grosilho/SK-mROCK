#ifndef SDEPROBLEMS_H
#define	SDEPROBLEMS_H

#include "DSde.h"
#include "OdeProblems.h"

// test equation
class DSdeDahlquistTestProblem: public DSde, public DahlquistTestProblem
{
public:
    DSdeDahlquistTestProblem();
    virtual ~DSdeDahlquistTestProblem(){};
 
    void g(Real t, Vector& x, Vector& G);
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
        
    Real phi(const Vector& X);
    
protected:
    Real sigma;
};

// Sde One dimensional test 1 from SROCK2 paper
class DSdeScalarNonStiffNonLinearTest: public DSde, public ScalarNonStiffNonLinearTest
{
public:
    DSdeScalarNonStiffNonLinearTest();
    virtual ~DSdeScalarNonStiffNonLinearTest(){};
 
    void g(Real t, Vector& x, Vector& G);
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
        
    Real phi(const Vector& X);
};


// Problem 2 in SROCK2 paper
class ManyDiffusionTerms: public DSde
{
public:
    ManyDiffusionTerms();
    virtual ~ManyDiffusionTerms();
    
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    void analytical_df(Real t, Vector& x, Matrix& dfx); 
    
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
    
    void rho(Real t, Vector& y, Real& eigmax);
        
    Real phi(const Vector& X);
        
protected:
    Vector a;
    Vector b;
};

/*
class SdeOneDimensionalTest3: public Sde
{
public:
    SdeOneDimensionalTest3(string output_file);
    virtual ~SdeOneDimensionalTest3(){};
 
    void init_solution();
    
    void f(Real t, Vector& x, Vector& fx);
    void g(Real t, Vector& x, Vector& gx);
        
    Real phi();
    Real Exact_phi();
    
    void rho(Real& eigmax);
protected:
    Real y0;
};

// Problem from Tempone paper.
// f is not differentiable so do not expect 2nd order convergence.
// On the other hand its very good for efficiency tests.
class SdeOneDimensionalTest4: public Sde
{
public:
    SdeOneDimensionalTest4(string output_file);
    virtual ~SdeOneDimensionalTest4(){};
 
    void init_solution();
    
    void f(Real t, Vector& x, Vector& fx);
    void g(Real t, Vector& x, Vector& gx);
        
    Real phi();
    Real Exact_phi();
    
    void rho(Real& eigmax);
protected:
    Real alpha;
};

// Problem from Tempone paper. 
class TwoDimensionalSinCos: public Sde
{
public:
    TwoDimensionalSinCos(string output_file);
    virtual ~TwoDimensionalSinCos(){};
 
    void init_solution();
    
    void f(Real t, Vector& x, Vector& fx);
    void g(Real t, Vector& x, Matrix& G);
    void g(Real t, Vector& x, Vector& G, int r);
        
    void rho(Real& eigmax);
    
    Real phi();
    Real Exact_phi();
};

class DuffingVanDerPol: public Sde
{
public:
    DuffingVanDerPol(string output_file);
    virtual ~DuffingVanDerPol(){};
 
    void init_solution();
    
    void f(Real t, Vector& x, Vector& fx);
    void g(Real t, Vector& x, Matrix& G);
    void g(Real t, Vector& x, Vector& G, int r);
            
    Real phi();
    Real Exact_phi();
    
    void rho(Real& eigmax);

protected:
    Real alpha;
    Real sigma;
};

class StochasticBrusselator: public Sde, public Brusselator
{
public:
    StochasticBrusselator(string output_file);
    virtual ~StochasticBrusselator(){};
 
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
        
    Real phi();
    Real Exact_phi();
    
protected:
    Real sigma;
};

class StochasticPopulationDynamics: public Sde, public PopulationDynamics
{
public:
    StochasticPopulationDynamics(string output_file);
    virtual ~StochasticPopulationDynamics(){};
     
    void g(Real t, Vector& x, Matrix& G);
    void g(Real t, Vector& x, Vector& G, int r);
        
    Real phi();
    Real Exact_phi();
        
protected:
    Real mu1;
    Real mu2;
};

class StochasticNeuronCable: public Sde, public NeuronCable
{
public:
    StochasticNeuronCable(string output_file);
    virtual ~StochasticNeuronCable(){};
     
    void g(Real t, Vector& x, Vector& G);
        
    Real phi();
    Real Exact_phi();
        
protected:
    Real sigma;
    Real V0;
};

*/
#endif	/* SDEPROBLEMS_H */

