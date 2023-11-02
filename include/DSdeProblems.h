#ifndef SDEPROBLEMS_H
#define SDEPROBLEMS_H

#include "DSde.h"

// test equation
class DahlquistTestProblem : public virtual DSde
{
public:
    DahlquistTestProblem();
    virtual ~DahlquistTestProblem(){};

    void set_initial_value(Vector &y0);
    void f(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmax);
    void AN_df(Real t, Vector &x, Matrix &dfx);

    void g(Real t, Vector &x, Vector &G);
    void g(Real t, Vector &x, Vector &G, int r);
    void g(Real t, Vector &x, Matrix &G);

    virtual Real phi(const Vector &X);

protected:
    static Real lambda;
    static Real mu;
    static Real xi;
};

// Problem of Section 5.1 in Abdulle, A., & Rosilho de Souza, G. (2022). Explicit stabilized multirate method for stiff stochastic differential equations. SIAM Journal on Scientific Computing, 44(4), A1859–A1883. https://doi.org/10.1137/21M1439018
class ScalarNonStiffNonLinearTest : public virtual DSde
{
public:
    ScalarNonStiffNonLinearTest();
    virtual ~ScalarNonStiffNonLinearTest(){};

    void set_initial_value(Vector &y0);
    void f(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmax);
    void AN_df(Real t, Vector &x, Matrix &dfx);

    void g(Real t, Vector &x, Vector &G);
    void g(Real t, Vector &x, Vector &G, int r);
    void g(Real t, Vector &x, Matrix &G);

    Real phi(const Vector &X);
};

// Problem of Section 5.2 in Abdulle, A., & Rosilho de Souza, G. (2022). Explicit stabilized multirate method for stiff stochastic differential equations. SIAM Journal on Scientific Computing, 44(4), A1859–A1883. https://doi.org/10.1137/21M1439018
class ManyDiffusionTerms : public virtual DSde
{
public:
    ManyDiffusionTerms();
    virtual ~ManyDiffusionTerms();

    void set_initial_value(Vector &y0);
    void f(Real t, Vector &x, Vector &fx);

    void g(Real t, Vector &x, Vector &G, int r);
    void g(Real t, Vector &x, Matrix &G);

    // void rho(Real t, Vector& y, Real& eigmax);

    Real phi(const Vector &X);

protected:
    Vector c;
    Vector R;
    Matrix Nu;
    Eigen::DiagonalMatrix<Real, 10> D;
};

class StochasticReactionDiffusion : public virtual DSde
{
public:
    StochasticReactionDiffusion();
    virtual ~StochasticReactionDiffusion();

    void set_initial_value(Vector &y0);

    void f(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmax);

    void g(Real t, Vector &x, Vector &G);

    Real phi(const Vector &X);

protected:
    Real nu;
    Real beta;
    Real sigma;
    Real V0;
};

class DiffusionRefinedMesh : public virtual DSde
{
public:
    DiffusionRefinedMesh();
    virtual ~DiffusionRefinedMesh();

    void set_initial_value(Vector &y0);

    void f(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmax);

    void g(Real t, Vector &x, Vector &G);

    Real phi(const Vector &X);

protected:
    Real nu;
    Real sigma;
    Real H1, H2, sqrtH1, sqrtH2;
    unsigned N1, N2;
};

class PopulationDynamics : public virtual DSde
{
public:
    PopulationDynamics();
    virtual ~PopulationDynamics();

    void set_initial_value(Vector &y0);

    void f(Real t, Vector &x, Vector &fx);

    void rho(Real t, Vector &y, Real &eigmax);

    void g(Real t, Vector &x, Vector &G, int r);
    void g(Real t, Vector &x, Matrix &G);

    Real phi(const Vector &X);

protected:
    Real l1, l2;
    Real mu1, mu2;
    Real alpha;
};

#endif /* SDEPROBLEMS_H */
