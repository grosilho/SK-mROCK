#ifndef MULTIRATESDEPROBLEMS_H
#define MULTIRATESDEPROBLEMS_H

#include "DSdeProblems.h"

class MultirateDahlquistTestProblem : public DahlquistTestProblem, public MultirateDSde
{
public:
    MultirateDahlquistTestProblem();
    virtual ~MultirateDahlquistTestProblem();

    void fF(Real t, Vector &x, Vector &fx);
    void fS(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS);
};

class MultirateScalarNonStiffNonLinearTest : public ScalarNonStiffNonLinearTest, public MultirateDSde
{
public:
    MultirateScalarNonStiffNonLinearTest();
    virtual ~MultirateScalarNonStiffNonLinearTest();

    void fF(Real t, Vector &x, Vector &fx);
    void fS(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS);
};

class MultirateManyDiffusionTerms : public ManyDiffusionTerms, public MultirateDSde
{
public:
    MultirateManyDiffusionTerms();
    virtual ~MultirateManyDiffusionTerms();

    void fF(Real t, Vector &x, Vector &fx);
    void fS(Real t, Vector &x, Vector &fx);

protected:
    Vector cF, cS;
    Vector FR, SR;
};

class MultirateStochasticReactionDiffusion : public StochasticReactionDiffusion, public MultirateDSde
{
public:
    MultirateStochasticReactionDiffusion();
    virtual ~MultirateStochasticReactionDiffusion();

    void fF(Real t, Vector &x, Vector &fx);
    void fS(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS);
};

class MultirateDiffusionRefinedMesh : public DiffusionRefinedMesh, public MultirateDSde
{
public:
    MultirateDiffusionRefinedMesh();
    virtual ~MultirateDiffusionRefinedMesh();

    void fF(Real t, Vector &x, Vector &fx);
    void fS(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS);
};

class MultiratePopulationDynamics : public PopulationDynamics, public MultirateDSde
{
public:
    MultiratePopulationDynamics();
    virtual ~MultiratePopulationDynamics();

    void fF(Real t, Vector &x, Vector &fx);
    void fS(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS);
};

#endif /* MULTIRATESDEPROBLEMS_H */