#ifndef ODEPROBLEMS_H
#define	ODEPROBLEMS_H

#include "Ode.h"
#include "lean_vtk.h"
#include "MyoKitIonicModel.h"
#include "Autodiff.h"

class DahlquistTestProblem: public virtual Ode
{
public:
    DahlquistTestProblem();
    virtual ~DahlquistTestProblem(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    void AN_df(Real t, Vector& x, Matrix& dfx); 
    void AD_df(Real t, Vector& x, Matrix& dfx); 
    static ADVector AD_f(Real t, ADVector& x);
        
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
    void AD_df(Real t, Vector& x, Matrix& dfx); 
    static ADVector AD_f(Real t, ADVector& x);
};

class NeuronCable: public virtual Ode
{
public:
    NeuronCable();
    virtual ~NeuronCable(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    void AD_df(Real t, Vector& x, SpMatrix& dfx); 
    static ADVector AD_f(Real t, ADVector& x);
    
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
    
    void AD_df(Real t, Vector& x, Matrix& dfx); 
    static ADVector AD_f(Real t, ADVector& x);
    
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
    
    void AD_df(Real t, Vector& x, SpMatrix& dfx); 
    static ADVector AD_f(Real t, ADVector& x);
        
protected:
    static Real A;
    static Real B;
    static Real alpha;
    static unsigned int Nu;
    static unsigned int Nv;
    static Real Hu;
    static Real Hv;
};


class Krogh10: public virtual Ode
{
public:
    Krogh10();
    virtual ~Krogh10(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    
    void AD_df(Real t, Vector& x, Matrix& dfx); 
    static ADVector AD_f(Real t, ADVector& x);
    
protected:
    static Real mu;
    static Real mus;
};

class PopulationDynamics: public virtual Ode
{
public:
    PopulationDynamics();
    virtual ~PopulationDynamics(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
        
protected:
    Real lambda1;
    Real lambda2;
    Real alpha;
};

class VanDerPol: public virtual Ode
{
public:
    VanDerPol();
    virtual ~VanDerPol(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
    
protected:
    Real eps;
};

class ODEIonicModel: public virtual Ode
{
public:
    ODEIonicModel();
    virtual ~ODEIonicModel(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& y, Vector& fy);
    
protected:
    unique_ptr<uBidomain::IonicModel> im;
    unsigned n_gating_vars;
    unsigned int npoints;
    Real CM;
    Real Iapp;
    Real interval;
    Real duration;
    Real first_stim;
    Real last_stim;
    
};



class DiffusionRefinedMesh: public virtual Ode
{
public:
    DiffusionRefinedMesh();
    virtual ~DiffusionRefinedMesh(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void rho(Real t, Vector& y, Real& eigmax);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
    
protected:
    Real nu;
    unsigned int N1,N2;
    Real H1,H2;
};

class InfectiousDiseaseTransmission: public virtual Ode
{
public:
    InfectiousDiseaseTransmission();
    virtual ~InfectiousDiseaseTransmission(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
    
protected:
    Real p,m,g,r,e,o;
    Real ID50, a;
};

// Taken from Hairer-Wanner II Ch. IV.1
class RobertsonChemicalSystem: public virtual Ode
{
public:
    RobertsonChemicalSystem();
    virtual ~RobertsonChemicalSystem(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
    void AD_df(Real t, Vector& x, Matrix& dfx); 
    static ADVector AD_f(Real t, ADVector& z);
    
protected:
    static Real k1,k2,k3;
};


class Oregonator: public virtual Ode
{
public:
    Oregonator();
    virtual ~Oregonator(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
    
protected:
    Real k1,k2,k3,k4,k5,ff,A,B;
};

class CUSP: public virtual Ode
{
public:
    CUSP();
    virtual ~CUSP(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
    
//protected:
    Real v(Real y);
    
protected:
    Real eps,sigma;
};

// Defines two problems taken from 
//  Yang et al, Spatial Resonances and Superposition Patterns in a Reaction-Diffusion Model with Interacting Turing Modes, 2002
//  Yang et al, Oscillatory Turing Patterns in Reaction-Diffusion Systems with Two Coupled Layers, 2003
// choose m=4 or m=5 in the constructor to solve one or the other problem
class mReactionDiffusion2DEquations: public virtual Ode
{
public:
    mReactionDiffusion2DEquations();
    virtual ~mReactionDiffusion2DEquations(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& y, Vector& fy);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
        
protected:
    Real inter_diff_and_reac(Vector& y, unsigned int k, unsigned int i, unsigned int j);
    
protected:
    unsigned int N, m, Nsq;
    vector<Real> D;
    Real alpha, beta, A, B;
    Real h;
};

class RadiationDiffusion: public virtual Ode
{
public:
    RadiationDiffusion();
    virtual ~RadiationDiffusion(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
    
protected:
    Real D1(Real x1, Real x2, Real E, Real gradE, Real T);
    Real D2(Real T);
    Real sigma(Real x1, Real x2, Real T);
    
    unsigned int Nsq, N;
    Real E0,T0;
    Real k;
    Real h;
};

class IntegroDifferentialEquation: public virtual Ode
{
public:
    IntegroDifferentialEquation();
    virtual ~IntegroDifferentialEquation(){};
 
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);
    void AN_df(Real t, Vector& x, Matrix& dfx);
    
//    static ADVector AD_f(Real t, ADVector& x, ADVector& fx);
    
protected:
    unsigned int N;
    Real sigma, h;
};

class MonoDomain: public virtual Ode
{
public:
    MonoDomain();
    virtual ~MonoDomain(){};
 
    virtual void set_initial_value(Vector& y0);
    virtual void f(Real t, Vector& x, Vector& fx);
        
protected:
    virtual Vector diff(const Vector& y);
    
    void init_mesh_data();
    void write_solution(const int nout, const string solname, 
                        const Real t, const Vector& y);

    unique_ptr<uBidomain::IonicModel> im;
    unsigned int N_gating_var;
    
    vector<Real> points;
    vector<int> elements;
    
    unsigned int Nx;
    unsigned int Ny;
    Real Lx;
    Real Ly;
    Real hx;
    Real hy;
    Real D_l;
    Real D_t;
    Real Cm;
    Real beta;
    Real scale_Iion;
    
    Vector diff_vect;
    Vector stim_vect;
};


#endif	/* ODEPROBLEMS_H */