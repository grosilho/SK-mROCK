#include "ClassicalOdeRungeKuttaIntegrators.h"
// #include "MatrixReplacement.h"

ExplicitEuler::ExplicitEuler(Parameters *param_, Ode *ode_)
    : OdeRungeKuttaIntegrator(param_, ode_)
{
    Astable = false;
    s = 1;

    reinit_integrator();
}

ExplicitEuler::~ExplicitEuler()
{
}

void ExplicitEuler::step(const Real t, const Real &h)
{
    Vector *&k1 = integr[0];

    ode->f(t, *yn, *k1);
    *ynpu = *yn + h * (*k1);

    n_f_eval++;
}

void ExplicitEuler::update_n_stages_and_h(Real &h)
{
    if (h > 0.9 * 2.0 / (this->eigmax))
    {
        h = 0.9 * 2.0 / (this->eigmax);
        this->last = false;
    }

    s_max = 1;
    s_avg += s;
}