#include "MultirateDSdeProblems.h"

MultirateDahlquistTestProblem::MultirateDahlquistTestProblem()
    : DahlquistTestProblem(), MultirateDSde()
{
}

MultirateDahlquistTestProblem::~MultirateDahlquistTestProblem()
{
}

void MultirateDahlquistTestProblem::fF(Real t, Vector &x, Vector &fx)
{
    fx(0) = lambda * x(0);
}

void MultirateDahlquistTestProblem::fS(Real t, Vector &x, Vector &fx)
{
    fx(0) = mu * x(0);
}

void MultirateDahlquistTestProblem::rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS)
{
    eigmaxF = abs(lambda);
    eigmaxS = abs(mu);
}

MultirateScalarNonStiffNonLinearTest::MultirateScalarNonStiffNonLinearTest()
    : ScalarNonStiffNonLinearTest(), MultirateDSde()
{
}

MultirateScalarNonStiffNonLinearTest::~MultirateScalarNonStiffNonLinearTest()
{
}

void MultirateScalarNonStiffNonLinearTest::fF(Real t, Vector &x, Vector &fx)
{
    fx(0) = 0.5 * sqrt(x(0) * x(0) + 1.);
}

void MultirateScalarNonStiffNonLinearTest::fS(Real t, Vector &x, Vector &fx)
{
    fx(0) = 0.25 * x(0);
}

void MultirateScalarNonStiffNonLinearTest::rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS)
{
    eigmaxF = abs(0.5 * y(0) / sqrt(y(0) * y(0) + 1));
    eigmaxS = 0.25;
}

MultirateManyDiffusionTerms::MultirateManyDiffusionTerms()
    : ManyDiffusionTerms(), MultirateDSde()
{
    FR.resize(Wsize);
    FR << 0, 0, 0, 0, 1, 0, 1, 0, 1, 0;
    SR = Eigen::VectorXd::Ones(Wsize) - FR;
    cF = FR.cwiseProduct(c);
    cS = SR.cwiseProduct(c);
}

MultirateManyDiffusionTerms::~MultirateManyDiffusionTerms()
{
}

void MultirateManyDiffusionTerms::fF(Real t, Vector &X, Vector &fx)
{
    R << X(1), X(2), X(2), X(4), X(0) * X(3), X(1), X(4) * (X(4) - 1) / 2, X(5), X(0) * X(5), X(6);
    fx = Nu * (cF.cwiseProduct(R));
}

void MultirateManyDiffusionTerms::fS(Real t, Vector &X, Vector &fx)
{
    R << X(1), X(2), X(2), X(4), X(0) * X(3), X(1), X(4) * (X(4) - 1) / 2, X(5), X(0) * X(5), X(6);
    fx = Nu * (cS.cwiseProduct(R));
}

MultirateStochasticReactionDiffusion::MultirateStochasticReactionDiffusion()
    : StochasticReactionDiffusion(), MultirateDSde()
{
}

MultirateStochasticReactionDiffusion::~MultirateStochasticReactionDiffusion()
{
}

void MultirateStochasticReactionDiffusion::fF(Real t, Vector &x, Vector &fx)
{
    fx(0) = nu * 2. * (x(1) - x(0)) * (neqn - 1.) * (neqn - 1.);
    fx(neqn - 1) = nu * 2. * (x(neqn - 2) - x(neqn - 1)) * (neqn - 1.) * (neqn - 1.);
    for (int i = 1; i < neqn - 1; i++)
    {
        fx(i) = nu * (x(i - 1) - 2. * x(i) + x(i + 1)) * (neqn - 1.) * (neqn - 1.);
    }
}

void MultirateStochasticReactionDiffusion::fS(Real t, Vector &x, Vector &fx)
{
    for (int i = 0; i < neqn; i++)
    {
        fx(i) = -beta * x(i);
        if (abs(i / (neqn - 1.) - 0.5) < 0.1)
            fx(i) += 5. * exp(1. - 1e-2 / (1e-2 - (i / (neqn - 1.) - 0 - 5) * (i / (neqn - 1.) - 0.5)));
    }
}

void MultirateStochasticReactionDiffusion::rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS)
{
    eigmaxF = nu * 4. * (neqn - 1.) * (neqn - 1.);
    eigmaxS = beta;
}

MultirateDiffusionRefinedMesh::MultirateDiffusionRefinedMesh()
    : DiffusionRefinedMesh(), MultirateDSde()
{
}

MultirateDiffusionRefinedMesh::~MultirateDiffusionRefinedMesh()
{
}

void MultirateDiffusionRefinedMesh::fF(Real t, Vector &x, Vector &fx)
{
    fx(0) = 0.;
    fx(neqn - 1) = nu * 2. * (x(neqn - 2) - x(neqn - 1)) / H2 / H2;
    fx(N1) = nu * (H2 * x(N1 - 1) - (H1 + H2) * x(N1) + H1 * x(N1 + 1)) / (0.5 * H1 * H2 * (H1 + H2));
    for (int i = 1; i < N1; i++)
        fx(i) = 0.;
    for (int i = N1 + 1; i < N1 + N2; i++)
        fx(i) = nu * (x(i - 1) - 2. * x(i) + x(i + 1)) / H2 / H2;
}

void MultirateDiffusionRefinedMesh::fS(Real t, Vector &x, Vector &fx)
{
    fx(0) = nu * 2. * (x(1) - x(0)) / H1 / H1; // Neumann bnd cond
    fx(neqn - 1) = 0.;
    fx(N1) = 0.;
    for (int i = 1; i < N1; i++)
        fx(i) = nu * (x(i - 1) - 2. * x(i) + x(i + 1)) / H1 / H1;
    for (int i = N1 + 1; i < N1 + N2; i++)
        fx(i) = 0.;
}

void MultirateDiffusionRefinedMesh::rho(Real t, Vector &y, Real &eigmaxF, Real &eigmaxS)
{
    eigmaxF = nu * 4. * max(1. / H1 / H1, 1. / H2 / H2);
    eigmaxS = nu * 4. * min(1. / H1 / H1, 1. / H2 / H2);
}

MultiratePopulationDynamics::MultiratePopulationDynamics()
    : PopulationDynamics(), MultirateDSde()
{
}

MultiratePopulationDynamics::~MultiratePopulationDynamics()
{
}

void MultiratePopulationDynamics::fF(Real t, Vector &X, Vector &fx)
{
    fx(0) = -l1 * X(0) * (1. - X(0));
    fx(1) = 0.;
}

void MultiratePopulationDynamics::fS(Real t, Vector &X, Vector &fx)
{
    fx(0) = alpha * (X(1) - 1.);
    fx(1) = -l2 * X(1) * (1. - X(1));
}

void MultiratePopulationDynamics::rho(Real t, Vector &X, Real &eigmaxF, Real &eigmaxS)
{
    eigmaxS = abs(l2 * (1. - 2 * X(1)));
    eigmaxF = abs(l1 * (1. - 2 * X(0)));
}