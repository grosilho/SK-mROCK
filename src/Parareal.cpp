#include "Parareal.h"
#include <sys/stat.h>
#include <chrono> 
using namespace std::chrono; 
#ifdef _OPENMP
#include <omp.h>
#endif

Parareal::Parareal(Parameters param_)
:param(param_)
{
    
}

Parareal::~Parareal()
{
}

void Parareal::integrate()
{
    #ifdef _OPENMP
    unsigned int num_omp_threads = omp_get_max_threads();
    #else
    unsigned int num_omp_threads = 1;
    #endif

    unsigned int iter = 24;
   
    param.verbose = false;
    
    Parameters outer_param = param;
    
    // init outer integrator and data
    outer_param.output_freq = -1;
    outer_param.rk_name = param.outer_rk_name;
    outer_param.dt = param.outer_dt;
    if(param.dtadap)
    {
        outer_param.rtol = 100*param.rtol;
        outer_param.atol = 100*param.atol;
    }
    Ode* ode;
    OdeRungeKuttaIntegrator* outer_integrator;
    outer_param.initOdeIntegration(outer_integrator,ode);
    
    // limit number of threads to number of coarse intervals T/coarse_dt.
    // i.e. we dont want to use so many threads that coarse_dt is reduced
    unsigned int n_coarse_intervals = ceil(ode->get_tend()/param.outer_dt);
    param.n_threads = min(param.n_threads, n_coarse_intervals);
    
    param.rk_name = param.inner_rk_name;
    param.dt = param.inner_dt;
    vector<Parameters> inner_params(param.n_threads,param);
    vector<OdeRungeKuttaIntegrator*> inner_integrators(param.n_threads);
    for(unsigned int n=0; n<param.n_threads;n++)
        inner_params[n].initOdeIntegration(inner_integrators[n],ode);
    string inner_out_file = param.output_path; //path+filename
    
    // init data
    outer_y.assign(param.n_threads+1, Vector::Zero(param.problem_size));
    inner_y.assign(param.n_threads+1, Vector::Zero(param.problem_size));
    outer_int_sol.assign(param.n_threads, Vector::Zero(param.problem_size));
    outer_t.resize(param.n_threads+1);
    Vector::Map(&outer_t[0], param.n_threads+1) = Vector::LinSpaced(param.n_threads+1,0.,ode->get_tend());
    
    ode->set_initial_value(outer_y[0]);
    ode->set_initial_value(inner_y[0]);
    
    Real elapsed_cpu_time = get_cpu_time();
    auto clock_time_start = high_resolution_clock::now();
    
    // first guess of coarse solution
    cout<<"Iteration 0: Compute first coarse guess"<<endl;
    for(unsigned int n=0;n<param.n_threads;n++)
    {
        cout<<"Call outer RK with t0 = "<<outer_t[n]<<", T = "<<outer_t[n+1]
            <<", out_dt = "<<param.outer_dt<<endl;
        
        outer_integrator->integrate(outer_y[n],outer_t[n],outer_t[n+1],
                                    outer_param.outer_dt,false);
        outer_y[n+1] = outer_integrator->solution();
        outer_int_sol[n] = outer_integrator->solution();
    }
    if(param.output_freq>=0)
        write_outer_y(ode,0);
    
    for(unsigned int k=1; k<=iter;k++)
    {
        cout<<"\nIteration "<<k<<": compute fine solutions"<<endl;
        //compute fine solutions
        #pragma omp parallel for
        for(unsigned int n=0;n<param.n_threads;n++)
        {
            cout<<"Call inner RK with t0 = "<<outer_t[n]<<", T = "<<outer_t[n+1]
                <<", in_dt = "<<inner_params[n].inner_dt<<endl;
            
            inner_params[n].output_path = inner_out_file + string("_in_sol_iter_")
                    +to_string(k)+string("_thread_")+to_string(n);
            
            inner_integrators[n]->integrate(outer_y[n],outer_t[n],outer_t[n+1],
                                            inner_params[n].inner_dt,true);
            inner_y[n+1] = inner_integrators[n]->solution();
        }
        
        cout<<"Iteration "<<k<<": correct outer solution"<<endl;
        // update the coarse solution
        for(unsigned int n=0;n<param.n_threads;n++)
        {
            cout<<"Call outer RK with t0 = "<<outer_t[n]<<", T = "<<outer_t[n+1]
                <<", out_dt = "<<param.outer_dt<<endl;

            outer_integrator->integrate(outer_y[n],outer_t[n],outer_t[n+1],
                                        outer_param.outer_dt,false);
            outer_y[n+1] = inner_y[n+1]
                           +outer_integrator->solution()-outer_int_sol[n];
            outer_int_sol[n] = outer_integrator->solution();
        }
        if(param.output_freq>=0)
            write_outer_y(ode,k);
    }
    
    elapsed_cpu_time = get_cpu_time();
    auto clock_time_stop = high_resolution_clock::now();
    auto elapsed_clock_time = 1e-6*duration_cast<microseconds>
                                  (clock_time_stop - clock_time_start).count();
    
    cout<<"\nCPU time: "<<elapsed_cpu_time<<endl;
    cout<<"Clock time: "<<elapsed_clock_time<<endl;
}

void Parareal::write_outer_y(Ode* ode, unsigned int iter)
{
    // First we write the integrating variable in a matlab readable file
    unsigned int N = ode->get_system_size();
    ofstream outfile;
        
    if(param.matlab_output)
    {
        outfile.open(param.output_path+string("_out_sol_iter_")
                    +to_string(iter)+string("_evolution.m"), ofstream::out);
       
        for(unsigned int nout=0;nout<param.n_threads+1;nout++)
        {
            outfile<<setprecision(16)<<"t("<<nout+1<<")="<<outer_t[nout]<<";"<<endl;
            outfile<<"y("<<nout+1<<",:)=[";
            for(int i=0;i<N-1;i++)
                outfile<<outer_y[nout](i)<<",";
            outfile<<outer_y[nout](N-1)<<"];"<<endl;
        }
        outfile.close();
    }
    
//    //Then we write it in a binary file
//    if(param->bin_output)
//    {
//        if(nout==1)
//            outfile.open(param->output_path+string("_evolution.bin"), ios::out | ios::binary);
//        else
//            outfile.open(param->output_path+string("_evolution.bin"), ios::out | ios::binary | ios::app);
//        outfile.write((char*)&t, sizeof(double));
//        outfile.write((char*)&(*y)(0), N*sizeof(double));
//        outfile.close();
//    }
//        
//    
//    //Finally we call the Ode class writing method, which implements problem specific output
//    ode->write_solution(nout, param->output_path, *y);
}

Real Parareal::get_cpu_time()
{
    return (Real)clock() / CLOCKS_PER_SEC;
}