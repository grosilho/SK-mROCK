#ifndef SDEPROBLEMS_H
#define	SDEPROBLEMS_H

#include "DSde.h"

// test equation
class DSdeDahlquistTestProblem: public virtual DSde
{
public:
    DSdeDahlquistTestProblem();
    virtual ~DSdeDahlquistTestProblem(){};

    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    void AN_df(Real t, Vector& x, Matrix& dfx); 
 
    void g(Real t, Vector& x, Vector& G);
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
        
    virtual Real phi(const Vector& X);
    
protected:
    static Real lambda;
    static Real mu;
    static Real xi;
};

// Sde One dimensional test 1 from SROCK2 paper
class DSdeScalarNonStiffNonLinearTest: public virtual DSde
{
public:
    DSdeScalarNonStiffNonLinearTest();
    virtual ~DSdeScalarNonStiffNonLinearTest(){};

    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    void AN_df(Real t, Vector& x, Matrix& dfx); 
 
    void g(Real t, Vector& x, Vector& G);
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
        
    Real phi(const Vector& X);
};


class ManyDiffusionTerms: public virtual DSde
{
public:
    ManyDiffusionTerms();
    virtual ~ManyDiffusionTerms();
    
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
    
    // void rho(Real t, Vector& y, Real& eigmax);
        
    Real phi(const Vector& X);
        
protected:
    Vector c;
    Vector R;
    Matrix Nu;
    Eigen::DiagonalMatrix<Real, 10> D;
};

class StochasticBrusselator: public virtual DSde
{
public:
    StochasticBrusselator();
    virtual ~StochasticBrusselator(){};
 
    void set_initial_value(Vector& y0);
    
    void f(Real t, Vector& x, Vector& fx);

    void g(Real t, Vector& x, Vector& G, int r);
    void g(Real t, Vector& x, Matrix& G);
        
    Real phi(const Vector& X);
    
protected:
    Real alpha;
    Real sigma;
};


class StochasticNeuronCable: public virtual DSde
{
public:
    StochasticNeuronCable();
    virtual ~StochasticNeuronCable();

    void set_initial_value(Vector& y0);
    
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
     
    void g(Real t, Vector& x, Vector& G);
        
    Real phi(const Vector& X);
        
protected:
    Real nu;
    Real beta;
    Real sigma;
    Real V0;
};


#endif	/* SDEPROBLEMS_H */

