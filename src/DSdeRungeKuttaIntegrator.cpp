#include "DSdeRungeKuttaIntegrator.h"

DSdeRungeKuttaIntegrator::DSdeRungeKuttaIntegrator(Parameters *param_, DSde *sde_)
    : OdeRungeKuttaIntegrator(param_, sde_), sde(sde_)
{
    needDoubleIntegral = false;

    if (sde->get_noise_type() != DIAGONAL)
        G.resize(sde->get_system_size(), sde->brownian_size());
}

DSdeRungeKuttaIntegrator::~DSdeRungeKuttaIntegrator()
{
}

void DSdeRungeKuttaIntegrator::reinit_statistics()
{
    OdeRungeKuttaIntegrator::reinit_statistics();

    n_g_eval = 0;
}

bool DSdeRungeKuttaIntegrator::integrate()
{
    HistoryStochasticIntegrals hsi(sde->brownian_size(),
                                   param->continuous, needDoubleIntegral, sde->get_noise_type(), param->seed);

    hsi.sample_integrals(sde->get_tend(), param->dt);

    return integrate(hsi);
}

bool DSdeRungeKuttaIntegrator::integrate(HistoryStochasticIntegrals &hsi)
{
    this->reinit_statistics();

    Vector *swap_ptr;
    bool fail_rho;
    nrho = 0;
    rho_outdated = true;
    sold = 0;

    sde->set_initial_value(*yn);
    Real h = param->dt;
    Real t = 0.;
    Real tend = sde->get_tend();
    Real told = t - 1.; // different from t

    elapsed_time = get_cpu_time();

    if (output_frequency > 0)
        output_solution(t, yn);

    for (unsigned int i = 0; i < hsi.length(); i++)
    {
        h = hsi[i].get_h();

        if (rho_outdated)
        {
            fail_rho = !update_rho(t);
            if (fail_rho)
                return false;
        }

        // The number of stages chosen depending on h and spectral radius.
        update_n_stages_and_h(h);

        I = &hsi[i];
        step(t, h);

        n_steps++; // update number of steps

        if (verbose)
            this->disp_step_info(t, h);

        nrho = nrho + 1;          // consecutive steps without computing the spectral radius
        nrho = (nrho) % rho_freq; // set to 0 every rho_freq steps
        rho_outdated = (nrho == 0);

        // ynpu becomes yn
        swap_ptr = yn;
        yn = ynpu;       // yn takes the new value
        ynpu = swap_ptr; // and ynpu will be free memory
        t = t + h;

        if ((last && output_frequency > 0) || (output_frequency > 0 && n_steps % output_frequency == 0))
            output_solution(t, yn);
    }

    elapsed_time = get_cpu_time() - elapsed_time;

    if (output_frequency >= 0)
        output_final_solution(yn);

    return true;
}

void DSdeRungeKuttaIntegrator::disp_step_info(Real &t, Real &h)
{
    cout << setprecision(4) << scientific;

    int stages = s + ((param->rk_name == "SROCK2") ? 2 : 0);
    std::string delta = u8"\u0394";

    cout << "Step t = " << setw(6) << setprecision(4) << t << ", " << delta << "t = " << setw(8) << setprecision(6) << h
         << ", s = " << setw(3) << stages << " and |y_n+1| = " << setw(7) << setprecision(4) << ynpu->lpNorm<Eigen::Infinity>() << ". " << endl;
}

void DSdeRungeKuttaIntegrator::convergence_test()
{
    // THIS CODE IS SERIAL
    // USE THE MONTECARLO CLASS IF YOU WANT TO RUN SIMULATIONS IN PARALLEL

    Vector refsol;
    Real tend = ode->get_tend();
    vector<Real> error(param->max_pow - param->min_pow + 1, 0.);
    vector<Real> rate(param->max_pow - param->min_pow + 1, 0.);

    for (unsigned int mc = 0; mc < param->MCiter; mc++)
    {
        //        cout<<"MC iter "<<mc+1<<endl;
        HistoryStochasticIntegrals hsi(sde->brownian_size(),
                                       param->continuous, needDoubleIntegral, sde->get_noise_type(), param->seed);
        hsi.sample_integrals(sde->get_tend(), (unsigned int)pow(2, param->max_pow + 2));
        this->integrate(hsi); // computing a ref solution but with the same scheme. Should use another more accurate method.
        refsol = this->solution();
        hsi.coarse();

        for (unsigned int k = param->max_pow; k >= param->min_pow; k--)
        {
            hsi.coarse();
            this->integrate(hsi);
            *yn -= refsol;
            error[k - param->min_pow] += yn->norm() / param->MCiter;
        }
    }

    cout << setprecision(6) << scientific;
    ofstream ofile(param->output_path + string("_conv_results.csv"), ofstream::out);
    ofile << "err, rate" << endl;

    for (unsigned int k = param->min_pow; k <= param->max_pow; k++)
    {
        if (k > param->min_pow)
            rate[k - param->min_pow] = log2(error[k - param->min_pow - 1] / error[k - param->min_pow]);
        ofile << error[k - param->min_pow] << ", " << rate[k - param->min_pow] << endl;
    }
    ofile.close();

    // printing results
    unsigned int prec = 4;
    unsigned int num_size = prec + 6;
    unsigned int cell_size = num_size + 1;
    unsigned int N = error.size();
    unsigned int table_width = 6 + (cell_size + 1) * N;

    unsigned int short_width = (table_width - 21);

    cout << setprecision(prec) << scientific;

    cout << setfill('-') << setw(short_width / 2) << "-"
         << " Convergence Results "
         << setw(short_width / 2 + short_width % 2) << "-" << endl;

    cout << "|Err |";
    for (unsigned int i = 0; i < N; i++)
        cout << setfill(' ') << setw(cell_size) << right << error[i] << "|";
    cout << endl;

    cout << setfill('-') << setw(table_width) << "-" << endl;
    cout << "|Rate|";
    for (unsigned int i = 0; i < N; i++)
        cout << setfill(' ') << setw(cell_size) << right << rate[i] << "|";
    cout << endl;
    cout << setfill('-') << setw(table_width) << "-" << endl;
}

void DSdeRungeKuttaIntegrator::print_integration_info()
{
    /**
     * Some statistics about the time integration.
     */

    cout << scientific;

    cout << "\n----------------   Integration Statistics   ----------------" << endl;

    if (Astable)
        cout << "The spectral radius has not been computed, " << param->rk_name << " is A-stable." << endl;
    else
    {
        cout << "Max estimation of the spectral radius: " << max_rho << endl;
        cout << "Min estimation of the spectral radius: " << min_rho << endl;
    }
    cout << "Number of f eval. for the spectr. radius = " << n_f_eval_rho << endl;
    cout << "Max number of stages used: " << s_max << endl;
    cout << "Mean number of stages used: " << ((Real)n_f_eval) / n_steps << endl;
    cout << "Number of f total evaluations = " << n_f_eval << endl;
    cout << "Number of g total evaluations = " << n_g_eval << endl;
    cout << "Time step used: " << h << endl;
    cout << "Steps: " << n_steps << endl;
    cout << "Number of steps: " << n_steps << endl;
    cout << "Elapsed time: " << elapsed_time << endl;
    cout << "------------------------------------------------------------\n"
         << endl;

    //    ofstream out(out_file+string("_statistics.txt"), ofstream::out);
    //    out<<setprecision(16);
    //    out<<elapsed_time<<", ";
    //    out<<n_f_eval<<", ";
    //    out<<max_s<<", ";
    //    out<<max_rho<<", ";
    //    out<<n_steps<<", ";
    //    out<<acc_steps<<", ";
    //    out<<rej_steps;
    //    out.close();
}

bool DSdeRungeKuttaIntegrator::need_double_integral()
{
    return needDoubleIntegral;
}

// -------------------------------------
MultirateDSdeRungeKuttaIntegrator::MultirateDSdeRungeKuttaIntegrator(Parameters *param_, MultirateDSde *msde_)
    : OdeRungeKuttaIntegrator(param_, msde_),
      DSdeRungeKuttaIntegrator(param_, msde_),
      MultirateOdeRungeKuttaIntegrator(param_, msde_),
      msde(msde_)
{
}

MultirateDSdeRungeKuttaIntegrator::~MultirateDSdeRungeKuttaIntegrator()
{
}

void MultirateDSdeRungeKuttaIntegrator::print_integration_info()
{
    /**
     * Some statistics about the time integration.
     */

    cout << scientific;

    string rho = u8"\u03C1";
    string delta = u8"\u0394";

    s_avg /= n_steps;
    m_avg /= n_steps;

    cout << "\n-------------------   Integration Info   -------------------" << endl;
    if (Astable)
        cout << "The spectral radius has not been computed, " << param->rk_name << " is A-stable." << endl;
    else
    {
        cout << "Max " << rho << "F: " << max_rho_F << endl;
        cout << "Min " << rho << "F: " << min_rho_F << endl;
        cout << "Max " << rho << "S: " << max_rho_S << endl;
        cout << "Min " << rho << "S: " << min_rho_S << endl;
        cout << "Number of fF eval. for " << rho << "F: " << n_fF_eval_rho << endl;
        cout << "Number of fS eval. for " << rho << "S: " << n_fS_eval_rho << endl;
    }

    cout << "Max s: " << s_max << endl;
    cout << "Mean s: " << s_avg << endl;
    cout << "Max m: " << m_max << endl;
    cout << "Mean m: " << m_avg << endl;
    cout << "fF evaluations = " << n_fF_eval << endl;
    cout << "fS evaluations = " << n_fS_eval << endl;
    cout << "g evaluations = " << n_g_eval << endl;
    cout << "Time step used: " << h << endl;
    cout << "Number of steps: " << n_steps << endl;
    cout << "Elapsed time: " << elapsed_time << endl;
    cout << "------------------------------------------------------------\n"
         << endl;
}

void MultirateDSdeRungeKuttaIntegrator::reinit_statistics()
{
    MultirateOdeRungeKuttaIntegrator::reinit_statistics();
    n_g_eval = 0;
}

void MultirateDSdeRungeKuttaIntegrator::disp_step_info(Real &t, Real &h)
{
    MultirateOdeRungeKuttaIntegrator::disp_step_info(t, h);
}