#include <cmath>
#include "DSdeProblems.h"

DahlquistTestProblem::DahlquistTestProblem()
    : Ode(), DSde()
{
    problem_name = "DahlquistTestProblem";

    tend = 0.1;

    neqn = 1;
    cte_rho = true;
    know_rho = true;
    analytical_df = true;
    dense_Jacobian = true;

    lambda = -100.;
    mu = -10;
    xi = 0.7 * sqrt(2. * abs(lambda + mu));

    Wsize = 1;
    noise = DIAGONAL;
}

Real DahlquistTestProblem::lambda;
Real DahlquistTestProblem::mu;
Real DahlquistTestProblem::xi;

void DahlquistTestProblem::set_initial_value(Vector &y0)
{
    y0(0) = 1.;
}

void DahlquistTestProblem::f(Real t, Vector &x, Vector &fx)
{
    fx(0) = (lambda + mu) * x(0);
}

void DahlquistTestProblem::AN_df(Real t, Vector &x, Matrix &dfx)
{
    dfx.resize(neqn, neqn);
    dfx(0, 0) = lambda + mu;
}

void DahlquistTestProblem::rho(Real t, Vector &y, Real &eigmax)
{
    eigmax = abs(lambda + mu);
}

Real DahlquistTestProblem::phi(const Vector &X)
{
    return sin(X(0));
}

void DahlquistTestProblem::g(Real t, Vector &x, Vector &G)
{
    G(0) = xi * x(0);
}

void DahlquistTestProblem::g(Real t, Vector &x, Vector &G, int r)
{
    G(0) = xi * x(0);
}

void DahlquistTestProblem::g(Real t, Vector &x, Matrix &G)
{
    G(0, 0) = xi * x(0);
}

// -----------------------------------------------------------------------------

ScalarNonStiffNonLinearTest::ScalarNonStiffNonLinearTest()
    : Ode(), DSde()
{
    problem_name = "ScalarNonStiffNonLinearTest";

    tend = 10.;

    neqn = 1;
    cte_rho = false;
    know_rho = true;
    analytical_df = true;
    dense_Jacobian = true;

    Wsize = 1;
    noise = DIAGONAL;
}

void ScalarNonStiffNonLinearTest::set_initial_value(Vector &y0)
{
    y0(0) = 0.;
}

void ScalarNonStiffNonLinearTest::f(Real t, Vector &x, Vector &fx)
{
    fx(0) = 0.25 * x(0) + 0.5 * sqrt(x(0) * x(0) + 1.);
}

void ScalarNonStiffNonLinearTest::AN_df(Real t, Vector &x, Matrix &dfx)
{
    dfx.resize(neqn, neqn);
    dfx(0, 0) = 0.25 + 0.5 * x(0) / sqrt(x(0) * x(0) + 1);
}

void ScalarNonStiffNonLinearTest::rho(Real t, Vector &y, Real &eigmax)
{
    eigmax = abs(0.25 + 0.5 * y(0) / sqrt(y(0) * y(0) + 1));
}

Real ScalarNonStiffNonLinearTest::phi(const Vector &X)
{
    return asinh(X(0));
}

void ScalarNonStiffNonLinearTest::g(Real t, Vector &x, Vector &G)
{
    G(0) = sqrt((x(0) * x(0) + 1.) / 2.);
}

void ScalarNonStiffNonLinearTest::g(Real t, Vector &x, Vector &G, int r)
{
    G(0) = sqrt((x(0) * x(0) + 1.) / 2.);
}

void ScalarNonStiffNonLinearTest::g(Real t, Vector &x, Matrix &G)
{
    G(0, 0) = sqrt((x(0) * x(0) + 1.) / 2.);
}

ManyDiffusionTerms::ManyDiffusionTerms()
    : Ode(), DSde()
{
    problem_name = "ManyDiffusionTerms";

    tend = 1.;

    neqn = 7;
    Wsize = 10;
    noise = GENERAL;

    cte_rho = false;
    know_rho = false;

    R.resize(Wsize);
    c.resize(Wsize);
    c << 0.0078, 0.043, 0.0039, 0.0007, 0.38, 3, 0.5, 5, 0.12, 9;

    Nu.resize(neqn, Wsize);
    Nu << 1., 0., 0., 0., -1., 1., 0., 0., -1., 1.,
        -1., 0., 0., 0., 1., -1., 0., 0., 0., 0.,
        1., 0., -1., 0., 0., 0., 0., 0., 0., 0.,
        1., 0., 0., 0., -1., 1., 0., 0., 0., 0.,
        0., 1., 0., -1., 0., 0., -2., 2., 0., 0.,
        0., 0., 0., 0., 0., 0, 1., -1., -1., 1.,
        0., 0., 0., 0., 0., 0., 0., 0., 1., -1.;
}

ManyDiffusionTerms::~ManyDiffusionTerms()
{
}

void ManyDiffusionTerms::set_initial_value(Vector &y0)
{
    y0.resize(neqn);
    y0 << 0., 1., 0., 30., 0., 0., 0.;
    y0 *= 5000.;
}

void ManyDiffusionTerms::f(Real t, Vector &X, Vector &fx)
{
    R << X(1), X(2), X(2), X(4), X(0) * X(3), X(1), X(4) * (X(4) - 1) / 2, X(5), X(0) * X(5), X(6);
    fx = Nu * (c.cwiseProduct(R));
}

void ManyDiffusionTerms::g(Real t, Vector &X, Vector &G, int r)
{
    R << X(1), X(2), X(2), X(4), X(0) * X(3), X(1), X(4) * (X(4) - 1) / 2, X(5), X(0) * X(5), X(6);
    G = Nu.col(r) * sqrt(abs(c(r) * R(r)));
}

void ManyDiffusionTerms::g(Real t, Vector &X, Matrix &G)
{
    R << X(1), X(2), X(2), X(4), X(0) * X(3), X(1), X(4) * (X(4) - 1) / 2, X(5), X(0) * X(5), X(6);
    D.diagonal() << c.cwiseProduct(R).cwiseAbs().cwiseSqrt();
    G = Nu * D;
}

Real ManyDiffusionTerms::phi(const Vector &X)
{
    return X(0) * X(0);
}

StochasticReactionDiffusion::StochasticReactionDiffusion()
    : Ode(), DSde()
{
    problem_name = "StochasticReactionDiffusion";

    tend = 1.;
    nu = 0.1;
    beta = 20.0;
    sigma = 2.;
    V0 = 10.;

    neqn = 128;
    cte_rho = true;
    know_rho = true;
    analytical_df = false;
    dense_Jacobian = false;

    Wsize = neqn;
    noise = DIAGONAL;
}

StochasticReactionDiffusion::~StochasticReactionDiffusion()
{
}

void StochasticReactionDiffusion::set_initial_value(Vector &y0)
{
    const Real pi = 4. * atan(1.);
    for (int j = 0; j < neqn; j++)
    {
        Real x = ((Real)j) / (neqn - 1.);
        y0(j) = -70. + 20. * cos(15. * pi * x) * (1. - x);
    }
}

void StochasticReactionDiffusion::f(Real t, Vector &x, Vector &fx)
{
    // Computing diffusion with Neumann bnd conditions
    fx(0) = nu * 2. * (x(1) - x(0)) * (neqn - 1.) * (neqn - 1.) - beta * x(0);
    fx(neqn - 1) = nu * 2. * (x(neqn - 2) - x(neqn - 1)) * (neqn - 1.) * (neqn - 1.) - beta * x(neqn - 1);
    for (int i = 1; i < neqn - 1; i++)
    {
        fx(i) = nu * (x(i - 1) - 2. * x(i) + x(i + 1)) * (neqn - 1.) * (neqn - 1.) - beta * x(i);
        if (abs(i / (neqn - 1.) - 0.5) < 0.1)
            fx(i) += 5. * exp(1. - 1e-2 / (1e-2 - (i / (neqn - 1.) - 0 - 5) * (i / (neqn - 1.) - 0.5)));
    }
}

void StochasticReactionDiffusion::rho(Real t, Vector &y, Real &eigmax)
{
    eigmax = nu * 4. * (neqn - 1) * (neqn - 1) + beta;
}

Real StochasticReactionDiffusion::phi(const Vector &X)
{
    return X.sum();
}

void StochasticReactionDiffusion::g(Real t, Vector &x, Vector &G)
{
    Real sqrtdx = 1. / sqrt(neqn - 1.);

    for (int i = 0; i < neqn; i++)
        G(i) = sigma * (x(i) + V0) / sqrtdx;
}

DiffusionRefinedMesh::DiffusionRefinedMesh()
    : Ode(), DSde()
{
    problem_name = "DiffusionRefinedMesh";

    tend = 1.;
    nu = 0.1;
    sigma = 1.;

    N1 = 100; // N1 points in [0,1]
    N2 = 200; // N2 points in ]1,2]

    if (N2 < N1)
    {
        cout << "Errore N2<N1" << endl;
        return;
    }

    neqn = N1 + N2 + 1;
    cte_rho = true;
    know_rho = true;
    analytical_df = false;
    dense_Jacobian = false;

    H1 = 1. / N1;
    H2 = 1. / N2;
    sqrtH1 = sqrt(H1);
    sqrtH2 = sqrt(H2);

    Wsize = neqn;
    noise = DIAGONAL;
}

DiffusionRefinedMesh::~DiffusionRefinedMesh()
{
}

void DiffusionRefinedMesh::set_initial_value(Vector &y0)
{
    const Real pi = 4. * atan(1.);
    Real x;
    for (int j = 0; j < N1; j++)
    {
        x = j * H1;
        y0(j) = sin(x * pi);
    }
    for (int j = N1; j <= N1 + N2; j++)
    {
        x = 1 + (j - N1) * H2;
        y0(j) = sin(x * pi);
    }
}

void DiffusionRefinedMesh::f(Real t, Vector &x, Vector &fx)
{
    fx(0) = nu * 2. * (x(1) - x(0)) / H1 / H1; // Neumann bnd cond
    fx(neqn - 1) = nu * 2. * (x(neqn - 2) - x(neqn - 1)) / H2 / H2;
    fx(N1) = nu * (H2 * x(N1 - 1) - (H1 + H2) * x(N1) + H1 * x(N1 + 1)) / (0.5 * H1 * H2 * (H1 + H2));
    for (int i = 1; i < N1; i++)
        fx(i) = nu * (x(i - 1) - 2. * x(i) + x(i + 1)) / H1 / H1;
    for (int i = N1 + 1; i < N1 + N2; i++)
        fx(i) = nu * (x(i - 1) - 2. * x(i) + x(i + 1)) / H2 / H2;
}

void DiffusionRefinedMesh::rho(Real t, Vector &y, Real &eigmax)
{
    eigmax = nu * 4. * max(1. / H1 / H1, 1. / H2 / H2);
}

Real DiffusionRefinedMesh::phi(const Vector &X)
{
    return X.sum();
}

void DiffusionRefinedMesh::g(Real t, Vector &x, Vector &G)
{
    for (int i = 0; i < N1; i++)
        G(i) = sigma * x(i) / sqrtH1;
    for (int i = N1 + 1; i < N1 + N2 + 1; i++)
        G(i) = sigma * x(i) / sqrtH2;
    G(N1) = sigma * x(N1) / sqrt((H1 + H2) / 2.);
}

PopulationDynamics::PopulationDynamics()
    : Ode(), DSde()
{
    problem_name = "PopulationDynamics";

    tend = 2.;

    neqn = 2;
    Wsize = 2;
    noise = GENERAL;

    cte_rho = false;
    know_rho = true;

    l1 = -50.;
    l2 = -1.;
    mu1 = 0.5 * sqrt(abs(l1));
    mu2 = 0.5 * sqrt(abs(l2));
    alpha = 2.;
}

PopulationDynamics::~PopulationDynamics()
{
}

void PopulationDynamics::set_initial_value(Vector &y0)
{
    y0.resize(neqn);
    y0 << 0.5, 0.5;
}

void PopulationDynamics::f(Real t, Vector &X, Vector &fx)
{
    fx(0) = -l1 * X(0) * (1. - X(0)) + alpha * (X(1) - 1.);
    fx(1) = -l2 * X(1) * (1. - X(1));
}

void PopulationDynamics::g(Real t, Vector &X, Vector &G, int r)
{
    if (r == 0)
    {
        G(0) = -mu1 * X(0) * (1. - X(0));
        G(1) = -mu2 * X(1) * (1. - X(1));
    }
    else if (r == 1)
    {
        G(0) = -mu2 * (1. - X(0));
        G(1) = 0.;
    }
}

void PopulationDynamics::g(Real t, Vector &X, Matrix &G)
{
    G(0, 0) = -mu1 * X(0) * (1. - X(0));
    G(0, 1) = -mu2 * (1. - X(0));
    G(1, 0) = -mu2 * X(1) * (1. - X(1));
    G(1, 1) = 0.;
}

Real PopulationDynamics::phi(const Vector &X)
{
    return X(0) * X(1);
}

void PopulationDynamics::rho(Real t, Vector &X, Real &eigmax)
{
    Real eigmaxS = abs(l2 * (1. - 2 * X(1)));
    Real eigmaxF = abs(l1 * (1. - 2 * X(0)));
    eigmax = max(eigmaxS, eigmaxF);
}
