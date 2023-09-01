#ifndef ODESTABILIZEDINTEGRATORS_H
#define	ODESTABILIZEDINTEGRATORS_H

#include "MainHeader.h"

#include "OdeRungeKuttaIntegrator.h"
#include "MultirateOdeRungeKuttaIntegrator.h"

class IERKC: public MultirateOdeRungeKuttaIntegrator
{
public:
    IERKC(Parameters* param_, MultirateOde* mode_);
    virtual ~IERKC();
    
protected:    
    void step(const Real t, const Real& h);
    
    void update_n_stages_and_h(Real& h);  
    void write_parameters(Real dt);
    
    virtual void disp_step_info(Real& t, Real& h, bool accepted);
    
protected:
    Vector* integr_add[4];
    Real damping, beta;
    unsigned int lin_solv_iter;
    Matrix J;
    Matrix I;
    SpMatrix spJ;
    SpMatrix spI;
    
    Eigen::BiCGSTAB<SpMatrix, Eigen::IncompleteLUT<Real>> solver; 
    Eigen::PartialPivLU<Matrix> dir_solver;
};

class mRKC: public MultirateOdeRungeKuttaIntegrator
{
public:
    mRKC(Parameters* param_, MultirateOde* mode_);
    virtual ~mRKC();
    
protected:    
    void step(const Real t, const Real& h);
    void f_eta(Real t, Vector& x, Vector& fx);
    
    void update_n_stages_and_h(Real& h);  
    void write_parameters(Real dt);
    
    virtual void disp_step_info(Real& t, Real& h, bool accepted);
    
protected:
    Vector* integr_add[4];
    Real damping, beta;
    Real eta;
};

class RKC1: public virtual OdeRungeKuttaIntegrator
{
public:
    RKC1(Parameters* param_, Ode* ode_);
    virtual ~RKC1();
    
protected:    
    virtual void step(const Real t, const Real& h);
    
    void reinit_integrator();
    
    void update_n_stages_and_h(Real& h);  
    void write_parameters(Real dt);
    void write_eigenvalues();
    
protected:
    Real damping, beta;
    unsigned int n_output_eigs;
};

class RKC2: public OdeRungeKuttaIntegrator
{
public:
    RKC2(Parameters* param_, Ode* ode_);
    virtual ~RKC2();
    
protected:    
    virtual void step(const Real t, const Real& h);
    
    void update_n_stages_and_h(Real& h);
    
    Real damping;
    Real beta;
};

class RKL1: public virtual OdeRungeKuttaIntegrator
{
public:
    RKL1(Parameters* param_, Ode* ode_);
    virtual ~RKL1();
    
protected:    
    virtual void step(const Real t, const Real& h);
    
    void update_n_stages_and_h(Real& h);  
    
protected:
    Real damping;
};

class RKL2: public OdeRungeKuttaIntegrator
{
public:
    RKL2(Parameters* param_, Ode* ode_);
    virtual ~RKL2();
    
protected:    
    virtual void step(const Real t, const Real& h);
    
    void update_n_stages_and_h(Real& h);
    
    Real damping;
};

class RKU1: public virtual OdeRungeKuttaIntegrator
{
public:
    RKU1(Parameters* param_, Ode* ode_);
    virtual ~RKU1();
    
protected:    
    virtual void step(const Real t, const Real& h);
    
    void update_n_stages_and_h(Real& h);  
    
protected:
    Real damping, beta;
};

class RKU2: public OdeRungeKuttaIntegrator
{
public:
    RKU2(Parameters* param_, Ode* ode_);
    virtual ~RKU2();
    
protected:    
    virtual void step(const Real t, const Real& h);
    
    void update_n_stages_and_h(Real& h);
    
    Real damping;
    Real beta;
};

class ROCK2: public OdeRungeKuttaIntegrator
{
public:
    ROCK2(Parameters* param_, Ode* ode_);
    virtual ~ROCK2();
    
protected:
    virtual void step(const Real t, const Real& h);

    virtual void update_n_stages_and_h(Real& h);
    void mdegr(int& mdeg, int mp[]);
    
 
protected:
    int mp[2];  ///<It is used in order to find the scheme's precomputed coefficients in the tables.
    
    static int ms[46];      ///<Array of coefficients.
    static Real fp1[46];    ///<Array of coefficients.
    static Real fp2[46];    ///<Array of coefficients.
    static Real recf[4476]; ///<Array of coefficients.
};

class DROCK2: public ROCK2
{
public:
    DROCK2(Parameters* param_, Ode* ode_);
    virtual ~DROCK2();
    
protected:
    virtual void step(const Real t, const Real& h);

    virtual void update_n_stages_and_h(Real& h);    
 
protected:
    static Real recf2[184]; ///<Array of coefficients.
    static Real recalph[46];///<Array of coefficients.
};


#endif	/* ODESTABILIZEDINTEGRATORS_H */