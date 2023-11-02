#include "Parameters.h"
#include "ClassicalOdeRungeKuttaIntegrators.h"
#include "StabilizedOdeRungeKuttaIntegrators.h"
#include "ClassicalDSdeRungeKuttaIntegrators.h"
#include "StabilizedDSdeRungeKuttaIntegrators.h"
#include "MultirateDSdeProblems.h"
#include <GetPot>
#include <filesystem>

Parameters::Parameters()
{
    problem_size = -1; // must be redefined later when initializing odes, sdes. If not an error will occur
    refsol_path = string("");
}

Parameters::~Parameters()
{
}

bool Parameters::initDSde(DSde *&sde)
{
    if (ntest == 1)
        sde = new MultirateDahlquistTestProblem();
    else if (ntest == 2)
        sde = new MultirateScalarNonStiffNonLinearTest();
    else if (ntest == 3)
        sde = new MultirateManyDiffusionTerms();
    else if (ntest == 4)
        sde = new MultirateStochasticReactionDiffusion();
    else if (ntest == 5)
        sde = new MultirateDiffusionRefinedMesh();
    else if (ntest == 6)
        sde = new MultiratePopulationDynamics();
    else
    {
        cout << "Problem " << ntest << " not known" << endl;
        return false;
    }

#pragma omp single
    {
        if (!std::filesystem::is_directory("../results/"))
            std::filesystem::create_directory("../results/");

        string output_dir = "../results/" + sde->get_problem_name();
        if (!std::filesystem::is_directory(output_dir))
            std::filesystem::create_directory(output_dir);

        output_path = "../results/" + sde->get_problem_name() + "/" + output_file;
        if (refsol_file.compare(string("")) != 0)
            refsol_path = "../results/" + sde->get_problem_name() + "/" + refsol_file;

        problem_size = sde->get_system_size();
    }

    return true;
}

bool Parameters::initDSdeTimeIntegrator(DSdeRungeKuttaIntegrator *&rk, DSde *sde)
{
    if (rk_name == "EM")
        rk = new EulerMaruyama(this, sde);
    else if (rk_name == "Platen")
        rk = new PlatenScheme(this, sde);
    else if (rk_name == "SKROCK")
        rk = new SKROCK(this, sde);
    else if (rk_name == "SKmROCK")
    {
        if (sde->is_multirate())
            rk = new SKmROCK(this, dynamic_cast<MultirateDSde *>(sde));
        else
        {
            cout << "ERROR: you cant use the multirate integrator " << rk_name << " with the non multirate problem " << sde->get_problem_name() << endl;
            return false;
        }
    }
    else
    {
        cout << "Integrator " << rk_name << " not known." << endl;
        return false;
    }

    return true;
}

bool Parameters::initDSdeIntegration(DSdeRungeKuttaIntegrator *&integr, DSde *&sde)
{
    if (!initDSde(sde))
        return false;

    return initDSdeTimeIntegrator(integr, sde);
}

bool Parameters::read_command_line(int argc, char **argv)
{
    GetPot cl(argc, argv); // command line parser

    ntest = cl.follow(ntest, 2, "-ntest", "-test");
    output_file = cl.follow(output_file.c_str(), 2, "-outputfile", "-ofile");
    refsol_file = cl.follow(refsol_file.c_str(), 3, "-refsol", "-ref", "-refsolfile");
    output_freq = cl.follow(output_freq, 2, "-outputfreq", "-ofreq");
    verbose = cl.follow(verbose, 2, "-verbose", "-verb");

    matlab_output = cl.follow(matlab_output, 3, "-matlab_output", "-matlab_out", "-matlab");
    bin_output = cl.follow(bin_output, 3, "-bin_output", "-bin_out", "-bin");

    rk_name = cl.follow(rk_name.c_str(), 2, "-solver", "-rk");
    dt = cl.follow(dt, 3, "-dt", "-tau", "-h");
    rho_freq = cl.follow(rho_freq, 3, "-rhofreq", "-rho_freq", "-rfreq");
    safe_add = cl.follow(safe_add, 2, "-safe_add", "-sa");

    if (cl.search("-convtest"))
        conv_test = true;
    else
        conv_test = false;
    max_pow = cl.follow(max_pow, 2, "-maxpow", "-max_pow");
    min_pow = cl.follow(min_pow, 2, "-minpow", "-min_pow");

    MCiter = cl.follow(MCiter, 3, "-iter", "-MCiter", "-mciter");
    continuous = cl.follow(continuous, "-contW");
    seed = cl.follow(seed, "-seed");
    if (cl.search("-pp"))
        post_proc = true;
    n_bins = cl.follow(n_bins, "-nbins");
    process_only = cl.follow(process_only, 3, "-process_only", "-processonly", "-po");

    if (MCiter > 1) // no verbose nor output if doing monte carlo
    {
        verbose = false;
        output_freq = -1;
    }

    if (conv_test)
    {
        verbose = false;
        output_freq = -1;
        rho_freq = 1;
    }

    if (cl.search("--help"))
    {
        cout << "This code implements the SK-ROCK and SK-mROCK explicit stabilized methods for stochastic differential equations.\n"
             << "Standard methods as Euler-Maruyama and Platen scheme are implemented for comparison.\n"
             << "Different problems are hardcoded in the executable. Look into src/DSdeProblems.cpp.\n"
             << "Run the code from the ./build folder as ./MultirateIntegrators OPTIONS, where OPTIONS is a combination of the below.\n"
                "Results are stored in the ./results folder.\n"
             << endl;
        cout << "The following options are available:\n"
             << "--- General options:\n"
             << "    -test       : a number in 1-6 specifying the problem that we want to solve. Default: 1\n"
             << "                  The list of problems is given below.\n"
             << "    -mciter     : Number of Monte Carlo iterations. Default: 1e3\n"
             << "    -ofile      : name of output file. Default: sol\n"
             << "    -rk         : The name of the numerical integrator to use. Default: SKROCK\n"
             << "                  The list of integrators is given below.\n"
             << "    -dt         : Time step size. Default: 1e-2\n"
             << "    -rfreq      : Frequency at which the spectral radii are re-estimated. Default: 5\n"
             << "    -safe_add   : Use some additional stages in explicit stabilized integrators. Default: 0\n"
             << "    -contW      : Generate Brownian motion from a continuous uniform distribution (true) \n"
             << "                  or a discrete distribution (false). Default: true\n"
             << "    -seed       : Fix the seed. Default: -1 \n"
             << "                  If -1 then a random seed is chosen.\n"
             << "    -verb       : Enables or disables verbosity. Default: true\n"
             << "--- Options for when doing many Monte Carlo iterations\n"
             << "    -ofreq      : Output frequency. Default: -1\n"
             << "                  - If 0 writes solution only at the end of simulation.\n"
             << "                    In general, used to generate a reference solution.\n"
             << "                  - If -1 never writes the solution.\n"
             << "                    In general, used to compute a solution and just compare it against a reference solution.\n"
             << "    -refsol     : name of the reference solution (if available). Default: \"\"\n"
             << "                  If available, errors are computed at the end of the simulation.\n"
             << "    -convtest   : If provided, performs a time convergence test, i.e. runs several simulations and checks errors. Default: false\n"
             << "                  If a reference solution is not provided, it is computed on the fly.\n"
             << "    -maxpow     : Minimal step size used for the convergence test is dt=tend/2^maxpow. Default: 6\n"
             << "    -minpow     : Maximal step size used for the convergence test is dt=tend/2^minpow. Default: 3\n"
             << "---- Options for when doing one Monte Carlo iteration (-mciter 1)\n"
             << "     Here are the options to write the solution on file, for later display with matlab scripts.\n"
             << "    -bin        : Writes solution in binary format. Default: false\n"
             << "    -matlab     : Writes solution in matlab format. Default: false\n"
             << "    -ofreq      : Output frequency. Default: -1\n"
             << "                  - If >0 writes solution every ofreq time steps,\n"
             << "--- List of problems:\n"
             << " This is the list of problems hardcoded in src/DSdeProblems.cpp. You choose them via the -test option."
             << " To change a parameter go to src/DSdeProblems.cpp, change it and recompile.\n"
             << " All problems are solvable with both a standard and a multirate solver.\n"
             << "     1 : The Dahlquist test equation.\n"
             << "     2 : A scalar non stiff non linear problem.\n"
             << "     3 : A problem with many diffusion terms.\n"
             << "     4 : A stochastic Reaction-Diffusion problem.\n"
             << "     5 : A heat equation with non uniform mesh.\n"
             << "     6 : A population dynamics problem.\n"
             << "--- List of numerical integrators:\n"
             << " Those marked with (m) are multirate integrators.\n"
             << "     EM     : Euler Maruyama.\n"
             << "     Platen : Platen Scheme.\n"
             << "     SKROCK : SK-ROCK method.\n"
             << "     SKmROCK: SK-mROCK multirate method.\n"
             << endl;

        return false;
    }

    return true;
}

void Parameters::print_info()
{
    cout << scientific;

    cout << "--------------- SDE Integration ---------------" << endl;
    cout << "Solver: " << rk_name << "." << endl;
    cout << "Step size: " << dt << "." << endl;
    cout << "Monte Carlo iterations: " << MCiter << "." << endl;
    cout << "Exact Brownian motion: " << (continuous ? "yes." : "no.") << endl;
    cout << "Seed: " << seed << "." << endl;
    if (conv_test)
    {
        cout << "Convergence test parameters: min_pow = " << min_pow << ", max_pow = " << max_pow << "." << endl;
    }
    else
    {
        cout << "Output file name: " << output_file << endl;
        cout << "Output frequency: " << output_freq << endl;
        cout << "Verbose: " << (verbose ? "yes." : "no.") << endl;
    }
    cout << "-----------------------------------------------" << endl;
}

void Parameters::print_info(DSde *sde)
{
    cout << scientific;

    cout << "----------------- SDE Problem -----------------" << endl;
    cout << "Problem name: " << sde->get_problem_name() << "." << endl;
    cout << "Problem size: " << sde->get_system_size() << "." << endl;
    cout << "Constant rho: " << (sde->is_rho_constant() ? "yes." : "no.") << endl;
    cout << "Known rho estimation: " << (sde->estimation_rho_known() ? "yes." : "no.") << endl;
    cout << "-----------------------------------------------" << endl;
}