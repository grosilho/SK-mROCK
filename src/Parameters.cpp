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
    problem_size = -1; //must be redefined later when initializing odes, sdes. If not an error will occur
    refsol_path = string("");
}        

Parameters::~Parameters()
{
}

bool Parameters::initDSde(DSde*& sde)
{
    if(ntest==1)
        sde = new MultirateDSdeDahlquistTestProblem();
    else if(ntest==2)
        sde = new MultirateDSdeScalarNonStiffNonLinearTest();
    else if(ntest==3)
        sde = new ManyDiffusionTerms();
//    else if(ntest==4)
//        sde = new SdeOneDimensionalTest4(output_file);
//    else if(ntest==5)
//        sde = new TwoDimensionalSinCos(output_file);
//    else if(ntest==6)
//        sde = new DuffingVanDerPol(output_file);
//    else if(ntest==7)
//        sde = new StochasticBrusselator(output_file);
//    else if(ntest==8)
//        sde = new StochasticPopulationDynamics(output_file);
//    else if(ntest==9)
//        sde = new StochasticNeuronCable(output_file);
    else
    { 
        cout<<"Problem "<<ntest<<" not known"<<endl;
        return false;
    }
    
    #pragma omp single
    {
        if (!std::filesystem::is_directory("../results/"))
        std::filesystem::create_directory("../results/");
    
        string output_dir = "../results/" +sde->get_problem_name();
        if (!std::filesystem::is_directory(output_dir))
            std::filesystem::create_directory(output_dir);
            
        output_path = "../results/" + sde->get_problem_name() + "/" + output_file;
        if(refsol_file.compare(string(""))!=0)
            refsol_path = "../results/" + sde->get_problem_name() + "/" + refsol_file;
        
        problem_size = sde->get_system_size();
    }
    
    return true;
}

bool Parameters::initDSdeTimeIntegrator(DSdeRungeKuttaIntegrator*& rk, DSde* sde)
{
    if(rk_name=="EM")
        rk = new EulerMaruyama(this,sde);
    else if(rk_name=="Platen")
        rk = new PlatenScheme(this,sde);
    else if(rk_name=="SKROCK")
        rk = new SKROCK(this,sde);
    else if(rk_name=="SKmROCK")
    {
        if(sde->is_multirate())
            rk = new SKmROCK(this,dynamic_cast<MultirateDSde*>(sde));
        else
        {
            cout<<"ERROR: you cant use the multirate integrator "<<rk_name<<" with the non multirate problem "<<sde->get_problem_name()<<endl;
            return false;
        }
    }
    else
    {
        cout<<"Integrator "<<rk_name<<" not known."<<endl;
        return false;
    }

    return true;
}

bool Parameters::initDSdeIntegration(DSdeRungeKuttaIntegrator*& integr, DSde*& sde)
{
    if(!initDSde(sde))
        return false;
    
    return initDSdeTimeIntegrator(integr,sde);
}


bool Parameters::read_command_line(int argc, char** argv)
{
    GetPot cl(argc, argv); //command line parser
    
    ntest = cl.follow(ntest,2,"-ntest","-test");
    output_file = cl.follow(output_file.c_str(),2,"-outputfile","-ofile");
    refsol_file = cl.follow(refsol_file.c_str(),3,"-refsol","-ref","-refsolfile");
    output_freq = cl.follow(output_freq,2,"-outputfreq","-ofreq");
    verbose = cl.follow(verbose,2,"-verbose","-verb");
    
    matlab_output = cl.follow(matlab_output,3,"-matlab_output","-matlab_out","-matlab");
    bin_output = cl.follow(bin_output,3,"-bin_output","-bin_out","-bin");
    specific_output = cl.follow(specific_output,3,"-spec_output","-spec_out","-spec");
    
    rk_name = cl.follow(rk_name.c_str(),2,"-solver","-rk");
    dt = cl.follow(dt,3, "-dt", "-tau","-h");
    rho_freq = cl.follow(rho_freq,3, "-rhofreq", "-rho_freq", "-rfreq");
    
    if(cl.search("-convtest"))
        conv_test=true;
    else
        conv_test=false;
    max_pow = cl.follow(max_pow,2,"-maxpow","-max_pow");
    min_pow = cl.follow(min_pow,2,"-minpow","-min_pow");
    
    MCiter = cl.follow(MCiter,3,"-iter","-MCiter","-mciter");
    continuous = cl.follow(continuous,"-contW");
    seed = cl.follow(seed,"-seed");
    if(cl.search("-pp"))
        post_proc=true;
    n_bins = cl.follow(n_bins,"-nbins");
    process_only = cl.follow(process_only,3,"-process_only","-processonly","-po");
    
    if(MCiter>1)//no verbose nor output if doing monte carlo
    {
        verbose=false;
        output_freq = -1;
    }
    
    if(conv_test)
    {
        verbose=false;
        output_freq=-1;
        rho_freq = 1;
    }
    
    if(cl.search("-help"))
    {
        cout<<"This is a code for running ODE and SDE simulations. A list of problems\n"
            <<"is hardcoded in the executable. Look into OdeList.h and SdeList.h.\n"
            <<"Each problem has a number n depending on the order of apparition in the list."<<endl;
        cout<<"The following options are available:\n"
            <<"    -sde       : run an SDE simulation.\n"
            <<"    -ode       : run an ODE simulation.\n"
            <<"    -oec       : sets the ODE error controller. Values are 1,2,3 for\n"
            <<"                 I=Integral, PI=Proportional I, PPI= Predictive PI controllers\n"
            <<"                 respectively. Default is 3.\n"
            <<"    -cee c     : if c=1 the integrator computes the exact local error (with a lot of Eulers steps).\n"
            <<"    -ewd e     : if e=1 then the error controller writes data it collects into disk at each step.\n"
            <<"    -ntest n   : chooses the problem. It takes the nth problem\n"
            <<"                 from list SdeList.h or from OdeList.h.\n"
            <<"    -dt h      : time step when running in fixed time step mode or initial time\n"
            <<"                 step when running in adaptive time step mode.\n"
            <<"    -dtadap dta: if dta=1 enables time step adaptivity, else in fixed time step mode.\n"
            <<"    -contW W   : if W=1 uses a continuous brownian motion, if W=0 a discrete one.\n"
            <<"    -atol at   : sets the absolute tolerance, by default at=1e-2.\n"
            <<"    -rtol rt   : sets the relative tolerance. If not provided then rtol=atol.\n"
            <<"    -ofile ofi : name of output file. By default ofi=solution.\n"
            <<"    -iter it   : Number of Monte Carlo simulations.\n"
            <<"    -ofreq ofr : frequency of write to disk. If ofr=0 writes at end of simulation\n"
            <<"                 only. If ofr=-1 never writes solution. If it>1 we set ofr=-1.\n"
            <<"    -verb v    : enables or disables verbosity. If it>1 then v=0.\n"
            <<"    -seed s    : sets the seed for random numbers. By default is -1, which is a random seed.\n"
            <<"    -solver rk : name of RungeKutta solver to use. The available solvers are:\n"
            <<"                 -RKC2     Runge-Kutta-Chebychev order 2,\n"
            <<"                 -ROCK2    Runge-Kutta-Orthogonal-Chebychev order 2,\n"
            <<"                 -DROCK2   ROCK2 with increased damping,\n"
            <<"                 -SROCK2   Stochastic Weak Order 2 ROCK2,\n"
            <<"                 -MT       Weak order 2 Milstein-Talay.\n"<<endl;
        cout<<"Error estimators are available for ROCK2 and SROCK2. Not implemented in DROCK2 and MT."<<endl;
        cout<<"Files are written into folders with same name of the tests."<<endl;
        return false;
    }
    
    return true;
}

void Parameters::print_info()
{
    cout<<scientific;
    
    cout<<"--------------- SDE Integration ---------------"<<endl;        
    cout<<"Solver: "<<rk_name<<"."<<endl;
    cout<<"Step size: "<<dt<<"."<<endl;    
    cout<<"Monte Carlo iterations: "<<MCiter<<"."<<endl;
    cout<<"Exact Brownian motion: "<<(continuous ? "yes.":"no.")<<endl;
    cout<<"Seed: "<<seed<<"."<<endl;
    if(conv_test)
    {
        cout<<"Convergence test parameters: min_pow = "<<min_pow<<", max_pow = "<<max_pow<<"."<<endl;
    }
    else
    {
        cout<<"Output file name: "<<output_path<<endl;
        cout<<"Output frequency: "<<output_freq<<endl;
        cout<<"Verbose: "<<(verbose ? "yes.":"no.")<<endl;                      
    }
    cout<<"-----------------------------------------------"<<endl;
    
}

void Parameters::print_info(DSde* sde)
{
    cout<<scientific;
        
    cout<<"----------------- SDE Problem -----------------"<<endl;
    cout<<"Problem name: "<<sde->get_problem_name()<<"."<<endl;
    cout<<"Problem size: "<<sde->get_system_size()<<"."<<endl;
    cout<<"Constant rho: "<<(sde->is_rho_constant() ? "yes.":"no.")<<endl;
    cout<<"Known rho estimation: "<<(sde->estimation_rho_known() ? "yes.":"no.")<<endl;
    cout<<"-----------------------------------------------"<<endl;
    
}