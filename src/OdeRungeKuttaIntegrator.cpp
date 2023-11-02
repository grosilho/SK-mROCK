#include <limits>
#include <cstdlib>
#include <string>
#include <fstream>

#include "OdeRungeKuttaIntegrator.h"
#include "StabilizedOdeRungeKuttaIntegrators.h"

OdeRungeKuttaIntegrator::OdeRungeKuttaIntegrator(Parameters *param_, Ode *ode_)
    : TimeIntegrator(param_, ode_)
{
    yn = new Vector(ode->get_system_size());
    ynpu = new Vector(ode->get_system_size());
    eigenvector = new Vector(ode->get_system_size());
    *eigenvector = Vector::Random(ode->get_system_size());

    for (int i = 0; i < 8; i++)
        integr[i] = new Vector(ode->get_system_size());

    Astable = false; // default

    verbose = param->verbose;
    output_frequency = param->output_freq;
    internal_rho = !ode->estimation_rho_known();
    rho_freq = param->rho_freq;

    cout << fixed << setfill(' ');

    //    reinit_integrator();
}

OdeRungeKuttaIntegrator::~OdeRungeKuttaIntegrator()
{
    delete yn;
    delete ynpu;
    delete eigenvector;
    for (int i = 0; i < 8; i++)
        delete integr[i];
}

void OdeRungeKuttaIntegrator::reinit_integrator()
{
    this->reinit_statistics();

    rho_outdated = true;
    nrho = 0;
    sold = 0;

    ode->set_initial_value(*yn);

    t = 0.;
    tend = ode->get_tend();
    h = param->dt;
}

void OdeRungeKuttaIntegrator::reinit_statistics()
{
    eigmax = 0.;
    max_rho = numeric_limits<int>::min();
    min_rho = numeric_limits<int>::max();
    s_max = 0;
    s_avg = 0.;
    n_f_eval_rho = 0;
    n_f_eval = 0;
    n_steps = 0;
    n_output = 1;
}

bool OdeRungeKuttaIntegrator::update_rho(Real t)
{
    /**
     * A new spectral radius is computed. Either with the estimation given by
     * ODE::rho or with the RungeKuttaIntegrator::rho internal power method.
     */

    if (Astable || (ode->is_rho_constant() && n_steps > 0))
        return true;

    if (verbose)
        cout << "\n--------------   Spectral Radius Estimation   --------------" << endl;

    int iter = 0;
    // Computed externally by ODE::rho
    if (!internal_rho)
        ode->rho(t, *yn, eigmax);
    // Computed internally by this->rho
    else
    {
        unsigned int conv_rho = this->rho(t, iter);
        if (conv_rho == 0)
        {
            cout << "ERROR: convergence failure in spectral radius computation. " << endl;
            return false;
        }
        else if (conv_rho == 1)
            cout << "WARNING: augment number of iterations in spectral radius computation." << endl;
    }

    // recover statistics
    max_rho = max(max_rho, (int)eigmax + 1);
    min_rho = min(min_rho, (int)eigmax + 1);

    nrho = 0;

    if (verbose)
        cout << scientific << "Spectral radius estimation: " << eigmax << endl;
    if (internal_rho && verbose)
        cout << "Power method converged in " << iter << " iterations." << endl;

    if (verbose)
        cout << "------------------------------------------------------------\n"
             << endl;

    rho_outdated = false;

    return true;
}

void OdeRungeKuttaIntegrator::disp_step_info(Real &t, Real &h)
{
    string rho = u8"\u03C1";
    cout << setprecision(4) << scientific;

    int stages = s + ((param->rk_name == "ROCK2" || param->rk_name == "DROCK2" || param->rk_name == "SROCK2") ? 2 : 0);
    std::string delta = u8"\u0394";

    cout << scientific;

    cout << "Step t = " << setw(6) << setprecision(4) << t << ", " << delta << "t = " << setw(8) << setprecision(6) << h
         << ", s = " << setw(3) << stages << ", " << rho << " = " << eigmax
         << " and |y_n+1| = " << setw(7) << setprecision(4)
         << ynpu->lpNorm<Eigen::Infinity>() << ". ";

    cout << endl;
}

Vector OdeRungeKuttaIntegrator::solution()
{
    return *yn;
}

unsigned int OdeRungeKuttaIntegrator::rho(Real t, int &iter)
{
    /**
     *     rho computes eigmax, a close upper bound of the
     *     spectral radius of the Jacobian matrix using a
     *     power method (J.N. Franklin (matrix theory)).
     *     The algorithm used is a small change (initial vector
     *     and stopping criteria) of that of
     *     Sommeijer-Shampine-Verwer, implemented in RKC.
     */

    Real eigmaxo, sqrtu, znor, ynor, quot, dzyn, dfzfn;

    const int maxiter = 100;
    const Real safe = 1.05;
    const Real tol = 1e-2;

    sqrtu = sqrt(numeric_limits<Real>::epsilon());

    // ------ The initial vectors for the power method are yn --------
    //       and yn+c*f(v_n), where vn=f(yn) a perturbation of yn
    //       (if n_steps=0) or a perturbation of the last computed
    //       eigenvector (if n_steps!=0).

    Vector *&fn = integr[0];
    Vector *z = integr[1];
    Vector *swap_ptr;

    // eigenvector is initialized as a random vector in the constructor
    // and then it contains an approx the the leading eigenvector
    *z = *eigenvector;

    ode->f(t, *yn, *fn);

    //    cout<<"norm yn "<<yn->norm()<<endl;
    //    cout<<"norm fn "<<fn->norm()<<endl;

    // ------ Perturbation.--------
    ynor = yn->norm();
    znor = z->norm();
    //    ynor= norm(*yn);
    //    znor= norm(*z);

    int k;

    // Building the vector z so that the difference z-yn is small
    if (ynor != 0.0 && znor != 0.0)
    {
        dzyn = ynor * sqrtu;
        quot = dzyn / znor;
        (*z) *= quot;
        (*z) += *yn;
    }
    else if (ynor != 0.0)
    {
        dzyn = ynor * sqrtu;
        *z = *yn;
        (*z) *= 1. + sqrtu;
    }
    else if (znor != 0.0)
    {
        dzyn = sqrtu;
        quot = dzyn / znor;
        (*z) *= quot;
    }
    else
    {
        dzyn = sqrtu * sqrt(z->size());
        for (int i = 0; i < ode->get_system_size(); i++)
            (*z)(i) += sqrtu;
    }
    // here dzyn=||z-yn|| and z=yn+(small perturbation)
    // dzyn=||z-yn|| will be always true, even with new z in the loop
    // Start the power method for non linear operator rhs

    //    eigmax=0.0;
    for (iter = 1; iter <= maxiter; iter++)
    {
        ode->f(t, *z, *eigenvector);
        n_f_eval_rho++;

        (*eigenvector) -= *fn; // dz is the new perturbation, not normalized yet
        dfzfn = eigenvector->norm();

        eigmaxo = eigmax;
        eigmax = dfzfn / dzyn; // approximation of the Rayleigh quotient (not with dot product but just norms)
        eigmax = safe * eigmax;

        if (abs(eigmax - eigmaxo) <= eigmax * tol)
        {
            // The last perturbation is stored. It will very likely be a
            //  good starting point for the next rho call.
            *eigenvector = *z;
            (*eigenvector) -= *yn;
            break;
        }
        if (dfzfn != 0.0)
        {
            quot = dzyn / dfzfn;
            *z = *eigenvector;
            (*z) *= quot;
            (*z) += *yn; // z is built so that dzyn=||z-yn|| still true
        }
        else
            return 0;
    }

    if (iter == maxiter + 1)
        return 1;
    else
        return 2;
}

void OdeRungeKuttaIntegrator::print_integration_info()
{
    string rho = u8"\u03C1";
    string delta = u8"\u0394";

    s_avg /= n_steps;

    cout << scientific;

    cout << "\n-------------------   Integration Info   -------------------" << endl;

    if (Astable)
        cout << "The spectral radius has not been computed, " << param->rk_name << " is A-stable." << endl;
    else
    {
        cout << "Max " << rho << ": " << max_rho << endl;
        cout << "Min " << rho << ": " << min_rho << endl;
        cout << "Number of f eval. for " << rho << ": " << n_f_eval_rho << endl;
    }
    cout << "Max s: " << s_max << endl;
    cout << "Mean s: " << s_avg << endl;
    cout << "f evaluations = " << n_f_eval << endl;
    cout << "Number of steps: " << n_steps << endl;
    cout << "Elapsed time: " << elapsed_time << endl;
    cout << "------------------------------------------------------------\n"
         << endl;
}

unsigned int OdeRungeKuttaIntegrator::get_n_f_eval()
{
    return n_f_eval;
}

unsigned int OdeRungeKuttaIntegrator::get_n_f_eval_rho_or_Nit()
{
    return n_f_eval_rho;
}

unsigned int OdeRungeKuttaIntegrator::get_n_steps()
{
    return n_steps;
}
