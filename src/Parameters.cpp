#include "Parameters.h"
#include "ClassicalOdeRungeKuttaIntegrators.h"
#include "StabilizedOdeRungeKuttaIntegrators.h"
#include "MultirateOdeProblems.h"
#include "ClassicalDSdeRungeKuttaIntegrators.h"
#include "StabilizedDSdeRungeKuttaIntegrators.h"
#include "DSdeProblems.h"
#include "ODECellByCellModel.h"
#include <GetPot>


Parameters::Parameters()
{
    problem_size = -1; //must be redefined later when initializind odes, sdes. If not an error will occur
    refsol_path = string("");
}        

Parameters::~Parameters()
{
}

bool Parameters::initOde(Ode*& ode)
{
    if(ntest==1)
        ode = new MultirateDahlquistTestProblem();
    else if(ntest==2)
        ode = new ScalarNonStiffNonLinearTest();
    else if(ntest==3)
        ode = new NeuronCable();
    else if(ntest==4)
        ode = new Brusselator();
    else if(ntest==5)
        ode = new MultiratePDEBrusselator();
    else if(ntest==6)
        ode = new Krogh10(); 
    else if(ntest==7)
        ode = new PopulationDynamics();
    else if(ntest==8)
        ode = new VanDerPol();
    else if(ntest==9)
        ode = new ODEIonicModel();
    else if(ntest==10)
        ode = new MultirateDiffusionRefinedMesh();
    else if(ntest==11)
        ode = new MultirateInfectiousDiseaseTransmission();
    else if(ntest==12)
        ode = new MultirateRobertsonChemicalSystem();
    else if(ntest==13)
        ode = new Oregonator();
    else if(ntest==14)
        ode = new MultirateCUSP();
    else if(ntest==15)
        ode = new MultiratemReactionDiffusion2DEquations();
    else if(ntest==16)
        ode = new MultirateRadiationDiffusion();
    else if(ntest==17)
        ode = new MultirateIntegroDifferentialEquation();
    else if(ntest==18)
        ode = new MultirateMonoDomain();
//    else if(ntest==19)
//        ode = new MonoDomain_in_1D();
    else if(ntest==20)
        ode = new MultirateODECellByCellModel(3,*this);//5e-4
    else
    {
        cout<<"Problem not known"<<endl;
        return false;
    }
    
    output_path = "../results/" + ode->get_problem_name() + "/" + output_file;
    if(refsol_file.compare(string(""))!=0)
        refsol_path = "../results/" + ode->get_problem_name() + "/" + refsol_file;
    
    problem_size = ode->get_system_size();
        
    return true;
}

bool Parameters::initOdeTimeIntegrator(OdeRungeKuttaIntegrator*& rk, Ode* ode)
{
    if(rk_name=="EE")
        rk = new ExplicitEuler(this,ode);
    else if(rk_name=="IE")
        rk = new ImplicitEuler(this,ode);
    else if(rk_name=="EM")
        rk = new ExplicitMidpoint(this,ode);
    else if(rk_name=="IM")
        rk = new ImplicitMidpoint(this,ode);
    else if(rk_name=="RK4")
        rk = new RungeKutta4(this,ode);
    else if(rk_name=="RKC1")
        rk = new RKC1(this,ode);
    else if(rk_name=="RKC2")
        rk = new RKC2(this,ode);
    else if(rk_name=="ROCK2")
        rk = new ROCK2(this,ode);
    else if(rk_name=="DROCK2")
        rk = new DROCK2(this,ode);
    else if(rk_name=="RKL1")
        rk = new RKL1(this,ode);
    else if(rk_name=="RKL2")
        rk = new RKL2(this,ode);
    else if(rk_name=="RKU1")
        rk = new RKU1(this,ode);
    else if(rk_name=="RKU2")
        rk = new RKU2(this,ode);
    else if(rk_name=="mRKC")
    {
        if(ode->is_multirate())
            rk = new mRKC(this,dynamic_cast<MultirateOde*>(ode));
        else
        {
            cout<<"ERROR: you cant use the multirate integrator "<<rk_name<<" with the non multirate problem "<<ode->get_problem_name()<<endl;
            return false;
        }
    }
    else if(rk_name=="IERKC")
    {
        if(ode->is_multirate())
            rk = new IERKC(this,dynamic_cast<MultirateOde*>(ode));
        else
        {
            cout<<"ERROR: you cant use the multirate integrator "<<rk_name<<" with the non multirate problem "<<ode->get_problem_name()<<endl;
            return false;
        }
    }
    
//    else if(rk_name=="SROCK2")
//    {
//        Sde* sde = dynamic_cast<Sde*>(ode);
//        rk = new SROCK2(sde, verbose, dtadap, atol, rtol, output_freq);
//    }
//    else if(rk_name=="MT")
//    {
//        Sde* sde = dynamic_cast<Sde*>(ode);
//        rk = new MilsteinTalay(sde, verbose, dtadap, atol, rtol, output_freq);
//    }
    else
    {
        cout<<"Integrator "<<rk_name<<" not known."<<endl;
        return false;
    }
    
    return true;
}

bool Parameters::initOdeIntegration(OdeRungeKuttaIntegrator*& rk, Ode*& ode)
{
    if(!initOde(ode))
        return false;
    
    return initOdeTimeIntegrator(rk,ode);
}

bool Parameters::initDSde(DSde*& sde)
{
    if(ntest==1)
        sde = new DSdeDahlquistTestProblem();
    else if(ntest==2)
        sde = new DSdeScalarNonStiffNonLinearTest();
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
    output_path = "../../../../Results/Tests/" + sde->get_problem_name() + "/" + output_file;
    if(refsol_file.compare(string(""))!=0)
        refsol_path = "../../../../Results/Tests/" + sde->get_problem_name() + "/" + refsol_file;
    
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
    
    if(cl.search("-ode"))
        eq = ODE;
    else if(cl.search(2,"-dsde","-d_sde"))
        eq = D_SDE;
    else if(cl.search(2,"-jdsde","-jd_sde"))
        eq = JD_SDE;
    
    rk_name = cl.follow(rk_name.c_str(),2,"-solver","-rk");
    dt = cl.follow(dt,3, "-dt", "-tau","-h");
    rho_freq = cl.follow(rho_freq,3, "-rhofreq", "-rho_freq", "-rfreq");
    
    if(cl.search("-parareal"))
        parareal=true;
    else
        parareal=false;
    n_threads = cl.follow(n_threads,2,"-n_threads","-n_cores");
    outer_rk_name = cl.follow(outer_rk_name.c_str(),3,"-outer_rk","-out_rk","-coarse_rk");
    inner_rk_name = cl.follow(inner_rk_name.c_str(),3,"-inner_rk","-in_rk","-fine_rk");
    outer_dt = cl.follow(outer_dt,3, "-outer_dt", "-out_dt","-out_h");
    inner_dt = cl.follow(inner_dt,3, "-inner_dt", "-in_dt","-in_h");
    
    if(cl.search("-convtest"))
        conv_test=true;
    else
        conv_test=false;
    max_pow = cl.follow(max_pow,2,"-maxpow","-max_pow");
    min_pow = cl.follow(min_pow,2,"-minpow","-min_pow");
    
    dtadap = cl.follow(dtadap,3,"-dtadaptivity","-dtadap","-dta");
    rtol = cl.follow(rtol,2,"-rtol","-rt");
    if(!cl.search(2,"-atol","-at"))
        atol=rtol;
    atol = cl.follow(atol,2,"-atol","-at");
    if(!cl.search(2,"-rtol","-rt"))
        rtol=atol;
    ode_contr = static_cast<Controller>(cl.follow(static_cast<int>(ode_contr),"-oec"));
    err_write_data = cl.follow(err_write_data,"-ewd");
    
    MCiter = cl.follow(MCiter,3,"-iter","-MCiter","-mciter");
    continuous = cl.follow(continuous,"-contW");
    seed = cl.follow(seed,"-seed");
    if(cl.search("-pp"))
        post_proc=true;
    n_bins = cl.follow(n_bins,"-nbins");
    process_only = cl.follow(process_only,3,"-process_only","-processonly","-po");
    
    //Some problem specific parameters
    P20_dG = cl.follow(1e-3,"-dG");
    P20_cell_l = cl.follow(0.,"-cl");
    P20_cell_w = cl.follow(0.,"-cw");
    P20_nx = cl.follow(0,"-nx");
    P20_ny = cl.follow(0,"-ny");
    P20_Rl = cl.follow(0.,"-Rl");
    P20_Rt = cl.follow(0.,"-Rt");
    P20_si = cl.follow(0.,"-si");
    P20_se = cl.follow(0.,"-se");
    P20_tend = cl.follow(0.,"-tend");
    P20_Imax = cl.follow(0.,"-Imax");
    P20_vertamp = cl.follow(0.,"-va");
    P20_vertfreq = cl.follow(0,"-vf");
    P20_vsmoothwave = cl.follow(true,"-vs");
    P20_horamp = cl.follow(0.,"-ha");
    P20_horlen = cl.follow(0.,"-hl");
    P20_horpos = cl.follow(0.,"-hp");
    P20_horprob = cl.follow(0.,"-hprob");
    P20_horalt = cl.follow(-2,"-halt");
    P20_hsmoothwave = cl.follow(true,"-hs");
    P20_onlyperp = cl.follow(false,"-op");
    P20_nic = cl.follow(1,"-nic");
    
    if(eq==ODE)//no monte carlo
        MCiter = 1;
    if(MCiter>1)//no verbose nor output if doing monte carlo
    {
        verbose=false;
        output_freq = -1;
    }
    if(eq!=ODE || rk_name=="EE" || rk_name=="IE" || rk_name=="EM" || rk_name=="IM" || rk_name=="RK4")
    {
        dtadap = false;
    }
    
    if(conv_test)
    {
        verbose=false;
        dtadap=false;
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
    
    if(eq==ODE)
    {
        cout<<"--------------- ODE Integration ---------------"<<endl;        
        cout<<"Solver: "<<rk_name<<"."<<endl;
        cout<<"Step size: "<<dt<<"."<<endl;
        cout<<"Step size adaptivity: "<<(dtadap ? "yes.":"no.")<<endl;
        if(dtadap)
        {
            cout<<"Controller type: ";
            if(ode_contr==I)
                cout<<"I."<<endl;
            else if(ode_contr==PI)
                cout<<"PI."<<endl;
            else if(ode_contr==PPI)
                cout<<"PPI."<<endl;
            cout<<"Relative tolerance: "<<rtol<<"."<<endl;
            cout<<"Absolute tolerance: "<<atol<<"."<<endl;
            cout<<"Write controller data: "<<(err_write_data ? "yes.":"no.")<<endl;
        }
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
    else
    {
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
}

void Parameters::print_info(Ode* ode)
{
    cout<<scientific;
    
    cout<<"----------------- ODE Problem -----------------"<<endl;
    cout<<"Problem name: "<<ode->get_problem_name()<<"."<<endl;
    cout<<"Problem size: "<<ode->get_system_size()<<"."<<endl;
    cout<<"Constant rho: "<<(ode->is_rho_constant() ? "yes.":"no.")<<endl;
    cout<<"Known rho estimation: "<<(ode->estimation_rho_known() ? "yes.":"no.")<<endl;
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