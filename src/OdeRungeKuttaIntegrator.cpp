#include <limits>
#include <cstdlib>
#include <string>
#include <fstream>

#include "OdeRungeKuttaIntegrator.h"
#include "ClassicalOdeRungeKuttaIntegrators.h"

OdeRungeKuttaIntegrator::OdeRungeKuttaIntegrator(Parameters* param_, Ode* ode_)
:TimeIntegrator(param_,ode_), err_control(param_)
{
    yn = new Vector(ode->get_system_size());
    ynpu = new Vector(ode->get_system_size());
    eigenvector = new Vector(ode->get_system_size());
    *eigenvector = Vector::Random(ode->get_system_size());
    
    for(int i=0;i<5;i++)
        integr[i]=new Vector(ode->get_system_size());
    
    uround=numeric_limits<Real>::epsilon();
            
    Astable = false; //default
    
    verbose=param->verbose;
    dt_adaptivity = param->dtadap;
    r_tol = param->rtol;
    a_tol = param->atol;
    output_frequency= param->output_freq;
    internal_rho = !ode->estimation_rho_known();
    rho_freq = param->rho_freq;
    
    err_order = 0; // must be redefined in subclasses
                
    cout<<fixed<<setfill(' ');
    
//    reinit_integrator();
}

 OdeRungeKuttaIntegrator::~OdeRungeKuttaIntegrator()
{
    delete yn;
    delete ynpu;
    delete eigenvector;
    for(int i=0;i<5;i++)
        delete integr[i];
}

void OdeRungeKuttaIntegrator::reinit_integrator()
{    
    this->reinit_statistics();
    
    rho_outdated=true;
    nrho = 0;
    sold=0;

    ode->set_initial_value(*yn);

    t = 0.;
    tend = ode->get_tend();   
    h = param->dt;
    
    err_control.reinit();
    err_control.update_hn(h);
    err_control.set_errD_order(err_order); 
}
 
void OdeRungeKuttaIntegrator::reinit_statistics()
{     
    eigmax = 0.;
    max_rho=numeric_limits<int>::min();
    min_rho=numeric_limits<int>::max();
    s_max=0;
    s_avg = 0.;
    dt_max = 0.;
    dt_min = numeric_limits<Real>::max();
    n_f_eval_rho=0;
    n_f_eval=0;
    n_steps = 0;
    acc_steps = 0;
    rej_steps = 0;
    n_output = 1;
}

bool  OdeRungeKuttaIntegrator::integrate(const Vector& y0, Real t0, Real T, 
                                        Real dt, bool new_integration)
{
    if(new_integration)
    {
        this->reinit_statistics();
        err_control.reinit();
        err_control.update_hn(h);
        err_control.set_errD_order(err_order); 
        rho_outdated=true;
        nrho = 0;
        sold=0;
    }
    //else we are continuing from the results of a previous call to this function
    
    *yn = y0;
    t = t0;
    tend = T;   
    h = dt;
            
    return integrate();
}

bool  OdeRungeKuttaIntegrator::integrate()
{

    /**
     * Ode is a class describing the Cauchy problem. It has a starting and end 
     * time step, time step, an initial value and a rhs function returning the 
     * right hand side.
     * This function integrates ode with a starting time step h. If one_step is
     * true just one step is executed. Recalling advance with the same arguments
     * does the next step using a new adapted h (id dt_adaptivity is true). 
     * If one_step is false integration is performed until end.
     * The returning value of idid is 1 if we reached the end, 2 if a step
     * has been executed correctly without reaching the end.
     */
    
    if(n_steps==0 && output_frequency>0)
        output_solution(t,yn);
    
    Vector *swap_ptr;
    bool fail_rho;
    
    last=false;
    
    elapsed_time = get_cpu_time();
    
    for(;;)
    {
        //if t+h>tend then adjust h
        if(t+1.01*h>=tend) 
        {
            h=tend-t;
            err_control.update_hn(h);
            last=true; //it will we be last step (if accepted)
        }
        if(h<10.0*uround) //too small h
        {
            cout<<"Step size is too small."<<endl;
            return false;  
        }
        
        if(rho_outdated) 
        {
            fail_rho = !(this->update_rho(t));
            if(fail_rho)
                return false;
        }

        //The number of stages chosen depending on h and spectral radius.
        update_n_stages_and_h(h);
        
        //Computation of an integration step.
        //this function is implemented in derived classes
        step(t,h);
        // at this point yn contains solution at time t, 
        // and ynpu is solution at time t+h
        // errors are already computed and passed to err_control
        
        n_steps++; //update number of steps 
       
        err_control.get_new_h(t, hnew);//compute hnew                

        if(verbose)
            this->disp_step_info(t,h,err_control.isAccepted());

        //Accepted step or without time step adaptivity
        if(err_control.isAccepted() || !dt_adaptivity)
        {
            //do some statistics
            dt_max = max(dt_max,h);
            dt_min = min(dt_min,h);
            
            //procedure for accepted step
            accepted_step(t,h);
            
            //ynpu becomes yn
            swap_ptr = yn; 
            yn = ynpu;      //yn takes the new value 
            ynpu = swap_ptr; //and ynpu will be free memory
                        
            if((last && output_frequency>0) || (output_frequency>0 && acc_steps%output_frequency==0))
                output_solution(t,yn);
            
            if(last)
                break;
        }
        else
            rejected_step(t,h);
    }
      
    elapsed_time = get_cpu_time()-elapsed_time;
    
    
    if(last && output_frequency>=0)
        output_final_solution(yn);
    
    if(last && param->refsol_path.compare(string(""))!=0)
        compute_errors();
    
    return last;
}     

 bool OdeRungeKuttaIntegrator::update_rho(Real t)
{
    /**
     * A new spectral radius is computed. Either with the estimation given by
     * ODE::rho or with the RungeKuttaIntegrator::rho internal power method.
     */
    
    if(Astable || (ode->is_rho_constant() && n_steps>0))
        return true;
    
    if(verbose) cout<<"\n--------------   Spectral Radius Estimation   --------------"<<endl;
        
    int iter=0;
    //Computed externally by ODE::rho
    if(!internal_rho)
        ode->rho(t,*yn,eigmax);
    //Computed internally by this->rho
    else
    {
        unsigned int conv_rho = this->rho(t,iter);
        if(conv_rho==0)
        {
            cout<<"ERROR: convergence failure in spectral radius computation. "<<endl;
            return false;
        }
        else if(conv_rho==1)
            cout<<"WARNING: augment number of iterations in spectral radius computation."<<endl;
    }

    //recover statistics
    max_rho = max(max_rho,(int)eigmax+1);
    min_rho = min(min_rho,(int)eigmax+1);        

    nrho=0;

    if(verbose) cout<<scientific<<"Spectral radius estimation: "<<eigmax<<endl;
    if(internal_rho && verbose)
        cout<<"Power method converged in "<<iter<<" iterations."<<endl;

    if(verbose) cout<<"------------------------------------------------------------\n"<<endl;
    
    rho_outdated=false;
    
    return true;
}

 void OdeRungeKuttaIntegrator::accepted_step(Real& t, Real& h)
{
    /**
     * Called when the step is accepted. Does output, sets the new time step size
     * and updated some variables.
     */

    acc_steps++;
    t=t+h;
    if(dt_adaptivity)
        h=hnew;//next h
    nrho=nrho+1; //consecutive steps without computing the spectral radius
    nrho=(nrho)%rho_freq; //set to 0 every rho_freq steps
    rho_outdated = (nrho==0);
}

void OdeRungeKuttaIntegrator::rejected_step(Real& t, Real& h)
{
    /**
     * Called when the step is rejected. Does output, chooses the new time step
     * and updates some variables.
     */
     
    rej_steps++;
    last=false;     //the step is rejected, it cant be the last one
    if(dt_adaptivity)
        h=hnew;

    //the step has been rejected, so we recompute the spectral radius, unless
    //it has just been computed (nrho==0)
    if(nrho>0)
        rho_outdated=true;
}

void OdeRungeKuttaIntegrator::disp_step_info(Real& t, Real& h, bool accepted)
{
    string rho =u8"\u03C1";
    cout<<setprecision(4)<<scientific;
    
    int stages = s+((param->rk_name=="ROCK2"||param->rk_name=="DROCK2"||param->rk_name=="SROCK2") ? 2:0);
    std::string delta = u8"\u0394";

    cout << scientific;
    
    cout<<"Step t = "<<setw(6)<<setprecision(4)<<t<<", "<<delta<<"t = "<<setw(8)<<setprecision(6)<<h
    <<", s = "<<setw(3)<<stages<<", "<<rho<<" = "<<eigmax
    <<" and |y_n+1| = "<<setw(7)<<setprecision(4)
    <<ynpu->lpNorm<Eigen::Infinity>()<<". ";
    
    if(accepted)
        cout<<" Accepted "; 
    else
        cout<<" Rejected ";
    err_control.disp_info();
    cout<<endl;
}

Vector OdeRungeKuttaIntegrator::solution()
{
    return *yn;
}

void OdeRungeKuttaIntegrator::convergence_test()
{
    Vector refsol;
    Real tend = ode->get_tend();
    vector<Real> error(param->max_pow-param->min_pow+1,0.);
    vector<Real> rate(param->max_pow-param->min_pow+1,0.);

    if(param->refsol_path.compare(string(""))!=0)
    {
        cout<<"Loading the reference solution."<<endl;
        read_reference_solution(&refsol);
    }
    else
    {
        cout<<"Computing a reference solution on the fly..."<<endl;
        verbose=true;
        param->dt = tend/pow(2,param->max_pow+2);
        TimeIntegrator* refintegrator = new RungeKutta4(param,ode);
        refintegrator->integrate();
        refsol = refintegrator->solution();
        delete refintegrator;
        verbose=false;
    }
    
    cout<<setprecision(6)<<scientific;
    ofstream ofile(param->output_path+string("_conv_results.csv"), ofstream::out);
    ofile<<"err, rate"<<endl;
    
    for(unsigned int k=param->min_pow;k<=param->max_pow;k++)
    {
        cout<<"Solving for k = "<<k<<", nsteps = "<<pow(2,k)<<", ";
        param->dt = tend/pow(2,k);
        this->reinit_integrator();
        this->integrate();
        *yn -= refsol;
        error[k-param->min_pow] = yn->norm()/refsol.norm();
        
        if(k>param->min_pow)
            rate[k-param->min_pow] = log2(error[k-param->min_pow-1]/error[k-param->min_pow]);
            
        cout<<"err = "<<error[k-param->min_pow]<<", rate = "<<rate[k-param->min_pow]<<"."<<endl;
        ofile<<error[k-param->min_pow]<<", "<<rate[k-param->min_pow]<<endl;
    }
    ofile.close();
    
    
    unsigned int prec = 4;
    unsigned int num_size = prec+6;
    unsigned int cell_size = num_size+1;
    unsigned int N = error.size();
    unsigned int table_width = 6+(cell_size+1)*N;
    
    unsigned int short_width = (table_width-21);
    
    cout<<setprecision(prec)<<scientific;
    
    cout<<setfill('-')<<setw(short_width/2)<<"-"<<" Convergence Results "
        <<setw(short_width/2+short_width%2)<<"-"<<endl;

    cout<<"|Err |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<error[i]<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    cout<<"|Rate|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<rate[i]<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
}


unsigned int OdeRungeKuttaIntegrator::rho(Real t, int& iter)
{
    /**
    *     rho computes eigmax, a close upper bound of the
    *     spectral radius of the Jacobian matrix using a 
    *     power method (J.N. Franklin (matrix theory)). 
    *     The algorithm used is a small change (initial vector
    *     and stopping criteria) of that of
    *     Sommeijer-Shampine-Verwer, implemented in RKC.
    */        
    
    Real eigmaxo,sqrtu,znor,ynor,quot,dzyn,dfzfn;
    
    const int maxiter=100;
    const Real safe=1.05;
    const Real tol = 1e-2;
    
    sqrtu= sqrt(uround);

// ------ The initial vectors for the power method are yn --------
//       and yn+c*f(v_n), where vn=f(yn) a perturbation of yn 
//       (if n_steps=0) or a perturbation of the last computed
//       eigenvector (if n_steps!=0). 
    
    Vector*& fn = integr[0];
    Vector* z  = integr[1];
    Vector* swap_ptr;
    
    //eigenvector is initialized as a random vector in the constructor
    //and then it contains an approx the the leading eigenvector
    *z=*eigenvector;     
        
    ode->f(t,*yn,*fn);
    
//    cout<<"norm yn "<<yn->norm()<<endl;
//    cout<<"norm fn "<<fn->norm()<<endl;
    
// ------ Perturbation.--------
    ynor= yn->norm();
    znor= z->norm();
//    ynor= norm(*yn);
//    znor= norm(*z);
    
    int k;
    
    // Building the vector z so that the difference z-yn is small
    if(ynor!=0.0 && znor!=0.0)
    {
        dzyn=ynor*sqrtu;
        quot=dzyn/znor;
        (*z) *= quot;
        (*z) += *yn;
    }
    else if(ynor!=0.0)
    {
        dzyn=ynor*sqrtu;
        *z=*yn;
        (*z) *= 1.+sqrtu;
    }
    else if(znor!=0.0)
    {
        dzyn=sqrtu;
        quot=dzyn/znor;
        (*z) *= quot;
    }
    else
    {
        dzyn=sqrtu*sqrt(z->size());
        for(int i=0;i<ode->get_system_size();i++)
            (*z)(i) += sqrtu;
    }
    //here dzyn=||z-yn|| and z=yn+(small perturbation)
    //dzyn=||z-yn|| will be always true, even with new z in the loop
    //Start the power method for non linear operator rhs    
    
//    eigmax=0.0;
    for(iter=1;iter<=maxiter;iter++)
    {
        ode->f(t,*z,*eigenvector);
        n_f_eval_rho++;

        (*eigenvector) -= *fn; //dz is the new perturbation, not normalized yet
        dfzfn= eigenvector->norm();
                
        eigmaxo=eigmax;
        eigmax=dfzfn/dzyn; //approximation of the Rayleigh quotient (not with dot product but just norms)
        eigmax=safe*eigmax;            
                        
        if(abs(eigmax-eigmaxo)<= eigmax*tol)
        {
            //The last perturbation is stored. It will very likely be a
            // good starting point for the next rho call.
            *eigenvector=*z;
            (*eigenvector) -= *yn;    
            break;
        }
        if (dfzfn!=0.0)
        {
            quot=dzyn/dfzfn;
            *z=*eigenvector;
            (*z) *= quot;
            (*z) += *yn; //z is built so that dzyn=||z-yn|| still true
        }
        else
            return 0;
    }
    
    if(iter==maxiter+1)
        return 1;
    else
        return 2;
}

void OdeRungeKuttaIntegrator::print_integration_info()
{    
    string rho =u8"\u03C1";
    string delta = u8"\u0394";
    
    s_avg /= n_steps;
    
    cout<<scientific;
    
    cout<<"\n-------------------   Integration Info   -------------------"<<endl;
    
    if(Astable)
        cout<<"The spectral radius has not been computed, "<<param->rk_name<<" is A-stable."<<endl;
    else
    {
        cout<<"Max "<<rho<<": "<<max_rho<<endl;
        cout<<"Min "<<rho<<": "<<min_rho<<endl;
        cout<<"Number of f eval. for "<<rho<<": "<<n_f_eval_rho<<endl;
    }
    cout<<"Max s: "<<s_max<<endl;
    cout<<"Mean s: "<<s_avg<<endl;
    cout<<"f evaluations = "<<n_f_eval<<endl;
    cout<<"Maximal "<<delta<<"t used: "<<dt_max<<endl;
    cout<<"Number of steps: "<<n_steps<<endl;
    cout<<"Accepted steps: "<<acc_steps<<endl;
    cout<<"Rejected steps: "<<rej_steps<<endl; 
    cout<<"Elapsed time: "<<elapsed_time<<endl;
    cout<<"------------------------------------------------------------\n"<<endl;
    
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

