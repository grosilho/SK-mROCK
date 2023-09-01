#include "MainHeader.h"
#include "Parameters.h"
#include "Parareal.h"
#include "MonteCarlo.h"


int main(int argc, char** argv)
{
  
//  Default parameters
    Parameters params; // Class for storing the parameters/options and initializing the integrators
    params.ntest = 20;           //problem number
    params.refsol_path = "";    //name of file containing reference solution. If empty, in a convergence test a reference solution is computed on the fly.
    params.verbose = true;      //show step info, otherwise just initial and final info.
    params.eq = ODE;            //ODE or D_SDE (SDE driven by diffusion) or JD_SDE (SDE driven by jump diffusion).
    params.rk_name = "RKC1";     //name of the solver
    params.dt = 1e-2;            //step size, or starting step size in case dtadap=true.
    params.rho_freq = 10;        //estimation of the spectral radius every rho_freq steps
    
    params.output_file = "sol"; //output file name
    params.output_freq = -1;    //-1 means no outputs, 0 just at the end, otherwise every output_freq steps
    params.matlab_output = false; // generates the .m file or not
    params.bin_output = false;    //generates the .bin file or not
    params.specific_output = false;//calls a problem specific output function
    
//  Parameters for convergence experiments
    params.conv_test = false;   //if false, we run once the experiment with step size params.dt
    params.max_pow = 7;         //if true we run a convergence test with step sizes tend/2^k,
    params.min_pow = 6;         //with k=min_pow,...,max_pow.

//  Parameters for parallel time integration
    params.parareal = false;
    params.n_threads = 10;
    params.outer_rk_name = "RKL1";
    params.inner_rk_name = "RKL1";    
    params.outer_dt = 0.5;
    params.inner_dt = 0.01;
    
//  Specific for step size adaptivity (only for ODEs)
    params.dtadap = false;      //toggles step size adaptivity.
    params.rtol = 1e-2;         //relative tolerance
    params.atol = 1e-2;         //absolute tolerance
    params.ode_contr = PPI;     //type of error controller: I (integral) or PPI (predictive proportional integral)
    params.err_write_data = false; //write error control data, as estimated errors, step sizes, predicted errors
    
//  Specific for SDEs
    params.MCiter = 10;         //Monte Carlo iterations
    params.continuous = true;   //Brownian motion (true) or discrete r.v. (false) with same first moments.
    params.seed = -1;           //if -1 the it's based on clock.
    params.post_proc=false;     //postprocessing for increasing ergodic order
    params.n_bins = 10;         //number of bins when computing density distance area
    params.process_only = false;//if true, the MonteCarlo class do not produce new data but import the two sets of data output_file and refsol. Processes and compares them.

    
    //read input parameters and eventually resolve incompatible choices. Exit if major error.
    if(!params.read_command_line(argc, argv))
        return 0;
    
    params.print_info();
        
    if(params.eq==ODE)
    {
        Ode* ode=0;
        OdeRungeKuttaIntegrator* integrator=0;
        
        if(!params.initOdeIntegration(integrator,ode))
            return 0;

        params.print_info(ode);
        
        if(params.parareal)
        {
            Parareal parareal(params);
            parareal.integrate();
        }
        else if(params.conv_test)
            integrator->convergence_test();
        else
        {
            integrator->integrate();
            integrator->print_integration_info();
        }
        
        delete ode;
        delete integrator;
          
    }
    else if(params.eq==D_SDE)
    {
        if(params.MCiter==1 && !params.process_only) // compute one sample
        {
            DSde* sde=0;
            DSdeRungeKuttaIntegrator* integrator=0;

            if(!params.initDSdeIntegration(integrator,sde))
                return 0;
            
            params.print_info(sde);

            integrator->integrate();
            integrator->print_integration_info();
            
            delete sde;
            delete integrator;
        }
        else//compute several samples and does some stats or does a convergence test
        {
            MonteCarlo montecarlo_integrator(params);
            
            if(params.conv_test)
                montecarlo_integrator.convergence_test();
            else
            {
                if(params.process_only)
                {
                    DSde* sde=0;
                    params.initDSde(sde);//only for updating the output_path variable
                    delete sde;
                    montecarlo_integrator.import_process_data();
                }
                else
                    montecarlo_integrator.integrate();
            }
        }
    }
    
    return 0;
}
