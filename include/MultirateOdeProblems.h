#ifndef FASTSLOWODEPROBLEMS_H
#define FASTSLOWODEPROBLEMS_H

#include "OdeProblems.h"


class MultirateDahlquistTestProblem: public DahlquistTestProblem, 
                                     public MultirateOde
{
public:
    MultirateDahlquistTestProblem();
    virtual ~MultirateDahlquistTestProblem();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);    
    
    void rho(Real t, Vector& y, Real& eigmax_F, Real& eigmax_S);
};

class MultiratePDEBrusselator: public PDEBrusselator, 
                              public MultirateOde
{
public:
    MultiratePDEBrusselator();
    virtual ~MultiratePDEBrusselator();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);    
};

//class MultirateODEIonicModel: public ODEIonicModel, 
//                              public MultirateOde
//{
//public:
//    MultirateODEIonicModel();
//    virtual ~MultirateODEIonicModel();
// 
//    void fF(Real t, Vector& y, Vector& fy);
//    void fS(Real t, Vector& y, Vector& fy);
//        
//};

class MultirateDiffusionRefinedMesh: public DiffusionRefinedMesh, 
                                     public MultirateOde
{
public:
    MultirateDiffusionRefinedMesh();
    virtual ~MultirateDiffusionRefinedMesh();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);
    
    void rho(Real t, Vector& y, Real& eigmax_F, Real& eigmax_S);
    
};

class MultirateInfectiousDiseaseTransmission: public InfectiousDiseaseTransmission, 
                                              public MultirateOde
{
public:
    MultirateInfectiousDiseaseTransmission();
    virtual ~MultirateInfectiousDiseaseTransmission();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);    
};

class MultirateRobertsonChemicalSystem: public RobertsonChemicalSystem, 
                                       public MultirateOde
{
public:
    MultirateRobertsonChemicalSystem();
    virtual ~MultirateRobertsonChemicalSystem();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);    
};

class MultirateCUSP: public CUSP, 
                    public MultirateOde
{
public:
    MultirateCUSP();
    virtual ~MultirateCUSP();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);   
};

class MultiratemReactionDiffusion2DEquations: public mReactionDiffusion2DEquations, 
                              public MultirateOde
{
public:
    MultiratemReactionDiffusion2DEquations();
    virtual ~MultiratemReactionDiffusion2DEquations();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);    
    
protected:
    unsigned int k_maxdiff;
};

class MultirateRadiationDiffusion: public RadiationDiffusion, 
                    public MultirateOde
{
public:
    MultirateRadiationDiffusion();
    virtual ~MultirateRadiationDiffusion();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);   
};

class MultirateIntegroDifferentialEquation: public IntegroDifferentialEquation, 
                    public MultirateOde
{
public:
    MultirateIntegroDifferentialEquation();
    virtual ~MultirateIntegroDifferentialEquation();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);   
};

class MultirateMonoDomain: public MonoDomain, 
                           public MultirateOde
{
public:
    MultirateMonoDomain();
    virtual ~MultirateMonoDomain();
 
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);   
};

#endif /* FASTSLOWODEPROBLEMS_H */

