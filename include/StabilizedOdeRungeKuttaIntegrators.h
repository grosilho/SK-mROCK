#ifndef ODESTABILIZEDINTEGRATORS_H
#define ODESTABILIZEDINTEGRATORS_H

#include "MainHeader.h"

#include "OdeRungeKuttaIntegrator.h"
#include "MultirateOdeRungeKuttaIntegrator.h"

class mRKC : public virtual MultirateOdeRungeKuttaIntegrator
{
public:
    mRKC(Parameters *param_, MultirateOde *mode_);
    virtual ~mRKC();

protected:
    void step(const Real t, const Real &h);
    void f_eta(Real t, Vector &x, Vector &fx);

    void update_n_stages_and_h(Real &h);
    void write_parameters(Real dt);

    virtual void disp_step_info(Real &t, Real &h);

protected:
    Vector *integr_add[4];
    Real damping, beta;
    Real eta;
    unsigned safe_add;
};

class RKC1 : public virtual OdeRungeKuttaIntegrator
{
public:
    RKC1(Parameters *param_, Ode *ode_);
    virtual ~RKC1();

protected:
    virtual void step(const Real t, const Real &h);

    void reinit_integrator();

    void update_n_stages_and_h(Real &h);
    void write_parameters(Real dt);
    void write_eigenvalues();

protected:
    Real damping, beta;
    unsigned safe_add;
    unsigned int n_output_eigs;
};

#endif /* ODESTABILIZEDINTEGRATORS_H */