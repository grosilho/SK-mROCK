#include "MonteCarlo.h"
#include <sys/stat.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static inline void loadBar(int x, int n, double start, double end);

MonteCarlo::MonteCarlo(Parameters param_)
:param(param_)
{
    DSde* sde=0;
    DSdeRungeKuttaIntegrator* integrator=0;
    param.initDSdeIntegration(integrator,sde);
    delete integrator;
    delete sde;
}

MonteCarlo::~MonteCarlo()
{
}

void MonteCarlo::integrate()
{  
    X.assign(param.MCiter, Vector::Zero(param.problem_size));
    vector<bool> succeed(param.MCiter,false);
    vector<Real> dt(param.MCiter,0.);
    vector<Real> n_f_eval(param.MCiter,0.);
    vector<Real> rho_Nit(param.MCiter,0.);
    
    #ifdef _OPENMP
    unsigned int num_threads = omp_get_max_threads();
    #else
    unsigned int num_threads = 1;
    #endif

    elapsed_time = get_cpu_time()/num_threads;
    
    #pragma omp parallel 
    {
        DSde* sde=0;
        DSdeRungeKuttaIntegrator* integrator=0;
        param.initDSdeIntegration(integrator,sde);
            
        #pragma omp for schedule(dynamic, 10)
        for(unsigned int i=0; i<param.MCiter; i++)
        {
            if(integrator->integrate()) // at each call samples a new Brownian path, on a mesh of dize param.dt
            {
                X[i] = integrator->solution();
                dt[i] = sde->get_tend()/integrator->get_n_steps();
                n_f_eval[i] = integrator->get_n_f_eval();
                rho_Nit[i] = integrator->get_n_f_eval_rho_or_Nit();
                succeed[i] = true;
            }
            loadBar(i, param.MCiter, elapsed_time, get_cpu_time()/num_threads);
        }
        
        delete integrator;
        delete sde;
    }
 
    elapsed_time = get_cpu_time()/num_threads-elapsed_time;
    

    process_data(succeed,dt,n_f_eval,rho_Nit); // also X
    write_data(succeed);
    print_data();
    
    if(param.refsol_path!=string(""))
        compute_print_errors();
}

void MonteCarlo::convergence_test()
{
    Vector strong_error = Vector::Zero(param.max_pow-param.min_pow+1);
    Vector weak_error = Vector::Zero(param.max_pow-param.min_pow+1);
    Vector strong_rate = Vector::Zero(param.max_pow-param.min_pow+1);
    Vector weak_rate = Vector::Zero(param.max_pow-param.min_pow+1);
 
    vector<Vector> strong_error_MC(param.MCiter,Vector::Zero(param.max_pow-param.min_pow+1));
    vector<Vector> weak_error_MC(param.MCiter,Vector::Zero(param.max_pow-param.min_pow+1));
    
    
    #ifdef _OPENMP
    unsigned int num_threads = omp_get_max_threads();
    #else
    unsigned int num_threads = 1;
    #endif

    elapsed_time = get_cpu_time()/num_threads;
        
    #pragma omp parallel 
    {
        DSde* sde=0;
        DSdeRungeKuttaIntegrator* integrator=0;
        param.initDSdeIntegration(integrator,sde);
        
        Vector refsol;
        Real refphi;
        HistoryStochasticIntegrals hsi(sde->brownian_size(),param.continuous,
                             integrator->need_double_integral(),sde->get_noise_type(),param.seed);
        //schedule(dynamic, 10)
        #pragma omp for schedule(dynamic, 10)
        for(unsigned int mc=0;mc<param.MCiter;mc++)
        {
            hsi.sample_integrals(sde->get_tend(),(unsigned int)pow(2,param.max_pow+2));
            integrator->integrate(hsi); // computing a ref solution but with the same scheme. Should use another more accurate method.
            refsol = integrator->solution();
            hsi.coarse();
            refphi = sde->phi(refsol);
            
            for(unsigned int k=param.max_pow;k>=param.min_pow;k--)
            {
                hsi.coarse();
                integrator->integrate(hsi);
                strong_error_MC[mc][k-param.min_pow] = Vector(integrator->solution()-refsol).norm();
                weak_error_MC[mc][k-param.min_pow] = sde->phi(integrator->solution())-refphi;
            }
            loadBar(mc, param.MCiter, elapsed_time, get_cpu_time()/num_threads);
        }
        
        delete integrator;
        delete sde;
    }
    elapsed_time = get_cpu_time()/num_threads-elapsed_time;
    
    cout<<"Elapsed time: "<<elapsed_time<<endl;
    
    for(unsigned int mc=0;mc<param.MCiter;mc++)
    {
        strong_error += strong_error_MC[mc];
        weak_error += weak_error_MC[mc];
    }
    strong_error/= param.MCiter;
    weak_error = weak_error.array().abs()/param.MCiter;
    
    
    cout<<setprecision(6)<<scientific;
    ofstream ofile(param.output_path+string("_conv_results.csv"), ofstream::out);
    ofile<<"strong_err, strong_rate, weak_err, weak_rate"<<endl;
    
    for(unsigned int k=param.min_pow;k<=param.max_pow;k++)
    {
        if(k>param.min_pow)
        {
            strong_rate[k-param.min_pow] = log2(strong_error[k-param.min_pow-1]/strong_error[k-param.min_pow]);
            weak_rate[k-param.min_pow] = log2(weak_error[k-param.min_pow-1]/weak_error[k-param.min_pow]);
        }
        ofile<<strong_error[k-param.min_pow]<<", "<<strong_rate[k-param.min_pow]<<", "<<weak_error[k-param.min_pow]<<", "<<weak_rate[k-param.min_pow]<<endl;
    }
    ofile.close();
    
    // printing results 
    unsigned int prec = 4;
    unsigned int num_size = prec+6;
    unsigned int cell_size = num_size+1;
    unsigned int N = strong_error.size();
    unsigned int table_width = 8+(cell_size+1)*N;
    
    unsigned int short_width = (table_width-21);
    
    cout<<setprecision(prec)<<scientific;
    
    
    cout<<setfill('-')<<setw(short_width/2)<<"-"<<" Convergence Results "
        <<setw(short_width/2+short_width%2)<<"-"<<endl;

    cout<<"|S Err |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<strong_error[i]<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    cout<<"|S Rate|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<strong_rate[i]<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;   
    
    cout<<"|W Err |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<weak_error[i]<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    cout<<"|W Rate|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<weak_rate[i]<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;   
}

void MonteCarlo::process_data(vector<bool>& succeed, vector<Real>& dt,
                 vector<Real>& n_f_eval, vector<Real>& rho_Nit)
{    
    cout<<"Processing data..."<<endl<<flush;
    
    mean = Array::Zero(param.problem_size);
    std_dev = Array::Zero(param.problem_size);
    dt_avg = 0.;
    n_f_eval_avg=0;
    rho_Nit_avg=0.;
    succeed_MCiter = 0;
    for(unsigned int i=0;i<param.MCiter;i++)
        if(succeed[i])
        {
            mean += X[i].array();
            dt_avg += dt[i];
            n_f_eval_avg += n_f_eval[i];
            rho_Nit_avg += rho_Nit[i];
            succeed_MCiter++;
        }
    mean /= succeed_MCiter;
    dt_avg /= succeed_MCiter;
    n_f_eval_avg /= succeed_MCiter;
    rho_Nit_avg /= succeed_MCiter;
    
    for(unsigned int i=0;i<param.MCiter;i++)
        if(succeed[i])
            std_dev += (X[i].array()-mean)*(X[i].array()-mean); // operations are element-wise, are we are using array instead of vector
    std_dev /= (succeed_MCiter-1.); 
    std_dev = std_dev.array().sqrt();
       
    printf("\033[F\033[J");
}

void MonteCarlo::print_data()
{
    if(param.refsol_path!=string(""))
        cout<<endl<<"----------------- SOLUTION DATA -----------------"<<endl<<endl;
    
    print_solver_info(elapsed_time,param.MCiter,succeed_MCiter,param.rk_name,
                      dt_avg,n_f_eval_avg,rho_Nit_avg,param.post_proc);
     
    print_statistics(mean, std_dev);
}

void MonteCarlo::print_solver_info(Real time, unsigned int iter, unsigned int succ_iter, 
                 string solver_name, Real dt_avg, Real n_f_eval_avg, Real rho_Nit_avg, 
                 bool postproc)
{
    cout.unsetf(ios::fixed | ios::scientific);
    cout<<setfill(' ');
    cout<<setprecision(4)<<scientific;
    
    cout<<"-------------- Solver info -------------"<<endl;
    cout<<setw(20)<<left<<"Solver: "<<solver_name<<endl;
//    cout<<setw(20)<<left<<"Mean nsteps : "<<param.get_final_time()/tau_mean<<endl;
    cout<<setw(20)<<left<<"Avg dt : "<<dt_avg<<endl;
    cout<<setw(20)<<left<<"Avg f eval : "<<n_f_eval_avg<<endl;
    cout<<setw(21)<<left<<"Avg \u03C1/Nit : "<<rho_Nit_avg<<endl;
//    cout<<setw(20)<<left<<"Newton tol : "<<Ntol<<endl;
    cout<<setw(20)<<left<<"Postprocessing: "<<(postproc? "yes":"no")<<endl;
    cout<<setw(20)<<left<<"Time to solution: "<<time<<" sec"<<endl;
    cout<<setw(20)<<left<<"Monte Carlo iter: "<<iter<<endl;
    cout<<setw(20)<<left<<"Successful iter: "<<succ_iter<<endl;
    cout.unsetf(ios::fixed | ios::scientific);
    cout<<setw(20)<<left<<"Succeed ratio: "<<100.*succ_iter/iter<<" %"<<endl;
    cout<<"----------------------------------------"<<endl;
}

void MonteCarlo::print_statistics(const Array& avg, const Array& std)
{
    unsigned int prec = 4;
    unsigned int num_size = prec+6;
    unsigned int cell_size = num_size+1;
    unsigned int N = avg.size();
    unsigned int table_width = 6+(cell_size+1)*N;
    
    unsigned int short_width = (table_width-12);
    
    cout<<setprecision(prec)<<scientific;
    
    
    cout<<setfill('-')<<setw(short_width/2)<<"-"<<" Statistics "
        <<setw(short_width/2+short_width%2)<<"-"<<endl;

    cout<<"|Avg |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<avg(i)<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    cout<<"|Std |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<std(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
}

void MonteCarlo::write_data(vector<bool>& succeed)
{
    cout<<"Writing data..."<<endl<<flush;
    
    // write all samples in a binary file
    unsigned int N = X[0].size();
    
    string filename = param.output_path+string("_data.bin");
    remove(filename.c_str());
    ofstream ofile(filename, ios::out | ios::binary);
    ofile.seekp(0);

    for(unsigned int i=0;i<param.MCiter;i++)
        if(succeed[i])
            ofile.write((char*)&X[i](0), N*sizeof(double));    
    ofile.close();
    
    // write statistics in a csv file
    filename = param.output_path+string("_statistics.csv");
    remove(filename.c_str());
    ofile.open(filename, ios::out);
    
    ofile<<setprecision(16);
    ofile<<"mean, std";
    for(unsigned int i=0;i<mean.size();i++)
        ofile<<endl<<mean(i)<<", "<<std_dev(i);
    ofile.close();
    
    // same as before but in a binary file
    filename = param.output_path+string("_statistics.bin");
    remove(filename.c_str());
    ofile.open(filename, ios::out | ios::binary);
    ofile.seekp(0);
    
    for(unsigned int i=0;i<mean.size();i++)
    {
        ofile.write((char*)&mean(i),sizeof(double));
        ofile.write((char*)&std_dev(i),sizeof(double));
    }
    ofile.close();
    
    // now write solver info in a csv file
    filename = param.output_path+string("_MCinfo.csv");
    remove(filename.c_str());
    ofile.open(filename, ios::out);
    
    ofile<<setprecision(16);
    ofile<<"time, MCiter, succeed_MCiter, solver, dt_avg, n_f_eval_avg, rho_Nit_avg, postproc"<<endl;
    ofile<<elapsed_time<<", "<<param.MCiter<<", "<<succeed_MCiter
         <<", "<<param.rk_name<<", "<<dt_avg<<", "<<n_f_eval_avg<<", "<<rho_Nit_avg<<", "<<param.post_proc;
    
    ofile.close();
    
    // same as before but in binary file
    filename = param.output_path+string("_MCinfo.bin");
    remove(filename.c_str());
    ofile.open(filename, ios::out | ios::binary);
    ofile.seekp(0);
    
    unsigned int solv_string_size =param.rk_name.length()+1;
    char solver_name[solv_string_size];
    strcpy(solver_name,param.rk_name.c_str());
    
    ofile.write((char*)&elapsed_time,sizeof(double));
    ofile.write((char*)&param.MCiter,sizeof(unsigned int));
    ofile.write((char*)&succeed_MCiter,sizeof(unsigned int));
    ofile.write((char*)&solv_string_size,sizeof(unsigned int)); //needed when reading the same file
    ofile.write((char*)&solver_name,sizeof(char)*solv_string_size);
    ofile.write((char*)&dt_avg,sizeof(double));
    ofile.write((char*)&n_f_eval_avg,sizeof(double));
    ofile.write((char*)&rho_Nit_avg,sizeof(double));
    ofile.write((char*)&param.post_proc,sizeof(bool));
    
    ofile.close();
    
    printf("\033[F\033[J");
}

void MonteCarlo::compute_print_errors()
{
    cout<<"Importing reference solution data..."<<endl<<flush;
    
    if(param.refsol_path.compare(string(""))==0)
    {
        cout<<"ERROR: there is no reference solution."<<endl;
        return;
    }
    
    //import reference solution
    vector<Vector> refX;
    import_data(param.refsol_path,refX);
    
    printf("\033[F\033[J");
    cout<<"Reading reference solution solver info..."<<endl<<flush;
    
    // reading SOLVER INFO of ref solution
    Real elapsed_time_ref;
    unsigned int MCiter_ref, succeed_MCiter_ref;
    Real dt_avg_ref,n_f_eval_avg_ref,rho_Nit_avg_ref;
    bool postproc_ref;
    string solver_name_ref;
    
    read_solver_info(param.refsol_path,elapsed_time_ref,MCiter_ref,succeed_MCiter_ref,
            solver_name_ref,dt_avg_ref,n_f_eval_avg_ref,rho_Nit_avg_ref,postproc_ref);
        
    printf("\033[F\033[J");
    cout<<"Reading reference solution statistics..."<<endl<<flush;
    
    // reading statistics of ref sol
    Array mean_ref, std_dev_ref;
    read_statistics(param.refsol_path,mean_ref,std_dev_ref);
        
    printf("\033[F\033[J");
    cout<<"Computing density distance area..."<<endl<<flush;
    
    vector<unsigned int> n_bins;
    Array dda;
    compute_density_distance_area(n_bins,dda,X,refX);
        
    printf("\033[F\033[J");
    cout<<"Estimating self distances..."<<endl<<flush;
    
    Array sd1, sd2;
    estimate_self_distances(n_bins, X.size(), refX.size(), sd1, sd2);
    
    
    printf("\033[F\033[J");
    
    // Printing info of reference solution
    cout<<endl<<"------------ REFERENCE SOLUTION DATA ------------"<<endl<<endl;
    print_solver_info(elapsed_time_ref,MCiter_ref,succeed_MCiter_ref,
            solver_name_ref,dt_avg_ref,n_f_eval_avg_ref,rho_Nit_avg_ref,postproc_ref);
    print_statistics(mean_ref, std_dev_ref);
     
    
    // printing errors
    cout<<endl<<"-------------------- ERRORS --------------------"<<endl<<endl;
    
    // print error of mean and std
    unsigned int N = param.problem_size;
    Array mean_err = (mean-mean_ref)/mean_ref;
    Array std_dev_err = (std_dev-std_dev_ref)/std_dev_ref;
    mean_err = mean_err.abs();
    std_dev_err = std_dev_err.abs();
    unsigned int prec = 4;
    unsigned int num_size = prec+6;
    unsigned int cell_size = num_size+1;
    unsigned int table_width = 6+(cell_size+1)*N;
    unsigned int short_width = (table_width-12);
    
    cout<<setfill('-')<<setw(short_width/2-2)<<"-"<<" Relative Error "
        <<setw(short_width/2+short_width%2-2)<<"-"<<endl;

    cout<<"|Avg |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<mean_err(i)<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    cout<<"|Std |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<std_dev_err(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    string filename = param.output_path+string("_errors.csv");
    remove(filename.c_str());
    ofstream ofile(filename, ios::out);
    for(unsigned int i=0;i<mean_err.size();i++)
    {
        ofile<<mean_err(i)<<", "<<std_dev_err(i)<<endl;
    }
    ofile.close();
    
    // printing the nbins, self distances, density distance area
    cout<<endl;
    cout<<setfill('-')<<setw(short_width/2-6)<<"-"<<" Densisty distance area "
        <<setw(short_width/2+short_width%2-6)<<"-"<<endl;

    cout<<"|nbin|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<n_bins[i]<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    cout<<"|SDX |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<sd1(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    cout<<"|SDrX|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<sd2(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    cout<<"|DDA |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<dda(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
}

void MonteCarlo::compute_density_distance_area(vector<unsigned int>& n_bins, Array& dda,
        vector<Vector>& X1, vector<Vector>& X2)
{
    unsigned int N = param.problem_size;
    n_bins.resize(N,param.n_bins); // default n bins
    dda.resize(N);
    
    for(unsigned int j=0;j<N;j++)
    {
        Real min_val = X1[0](j);
        Real max_val = X1[0](j);
        
        #pragma omp parallel reduction(min:min_val) reduction(max:max_val)
        {
        #pragma omp for
        for(unsigned int i=0;i<X1.size();i++)
        {
            min_val = min(min_val,X1[i](j));
            max_val = max(max_val,X1[i](j));
        }
        #pragma omp for
        for(unsigned int i=0;i<X2.size();i++)
        {
            min_val = min(min_val,X2[i](j));
            max_val = max(max_val,X2[i](j));
        }
        }
        
        // adaptive bins selection
        min_val = floor(min_val);
        if(max_val != ceil(max_val))
            max_val = ceil(max_val);
        else
            max_val++;
        if((int)round(max_val-min_val)%2==1)
            max_val++;
        
        if(n_bins[j]==0)
        {
            n_bins[j] = max_val-min_val; //bins of size 1 integer
            //but many variables increase or decrease their value by two, so its better to have bins of size 2
            n_bins[j] = n_bins[j]/2;
        }// else use the default (not zero) value
        
        dda(j)=0.;
        Real nX1, nX2;
        Real bin_min, bin_max;
        for(unsigned int b=0;b<n_bins[j];b++)
        {
            nX1=0.;
            nX2=0.;
            bin_min = min_val + b*((max_val-min_val)/n_bins[j]);
            bin_max = min_val + (b+1.)*((max_val-min_val)/n_bins[j]);
            
            #pragma omp parallel reduction(+:nX1,nX2)
            {
            #pragma omp for
            for(unsigned int i=0;i<X1.size();i++)
                if(bin_min<=X1[i](j) && X1[i](j)<bin_max)
                    nX1++;
            #pragma omp for
            for(unsigned int i=0;i<X2.size();i++)
                if(bin_min<=X2[i](j) && X2[i](j)<bin_max)
                    nX2++;
            }

            dda(j) += abs(nX1/X1.size()-nX2/X2.size());
        }
    }
    
}

void MonteCarlo::estimate_self_distances(vector<unsigned int>& n_bins, unsigned int N1, unsigned int N2,
        Array& sd1, Array& sd2)
{
    unsigned int N = n_bins.size();
    sd1.resize(N);
    sd2.resize(N);
    for(unsigned int i=0;i<N;i++)
    {
        sd1(i) = sqrt(4.*n_bins[i]/3.141592/N1);
        sd2(i) = sqrt(4.*n_bins[i]/3.141592/N2);
    }
}

void MonteCarlo::import_process_data()
{  
    import_data(param.output_path,X);
    read_solver_info(param.output_path, elapsed_time, param.MCiter, succeed_MCiter, 
            param.rk_name,dt_avg,n_f_eval_avg,rho_Nit_avg,param.post_proc);
    read_statistics(param.output_path,mean,std_dev);
    
    print_data();
    
    compute_print_errors();
}

void MonteCarlo::import_data(string filename, vector<Vector>& Y)
{
    //import solution data
    filename = filename + string("_data.bin");

    struct stat st;
    stat(filename.c_str(), &st);
    unsigned int file_size = st.st_size;
    
    unsigned int N = param.problem_size;
    unsigned int samples = file_size/sizeof(double)/N;
    Y.resize(samples, Vector::Zero(N));
    
    ifstream ifile(filename, ios::in | ios::binary);
    for(unsigned int i=0;i<samples;i++)
            ifile.read((char*)&Y[i](0), N*sizeof(double));
    ifile.close();
}

void MonteCarlo::read_solver_info(string filename, Real& time, unsigned int& iter,
                 unsigned int& succeed_iter, string& solver_name_str,
                 Real& dt_avg, Real& n_f_eval_avg, Real& rho_Nit_avg, 
                 bool& post_proc)
{
    unsigned int solv_string_size;
    filename = filename+string("_MCinfo.bin");
    ifstream ifile(filename, ios::in | ios::binary);
    
    ifile.read((char*)&time,sizeof(double));
    ifile.read((char*)&iter,sizeof(unsigned int));
    ifile.read((char*)&succeed_iter,sizeof(unsigned int));
    ifile.read((char*)&solv_string_size,sizeof(unsigned int)); 
    char solver_name_char[solv_string_size];
    ifile.read((char*)&solver_name_char,sizeof(char)*solv_string_size);
    solver_name_str = solver_name_char;
    ifile.read((char*)&dt_avg,sizeof(double));
    ifile.read((char*)&n_f_eval_avg,sizeof(double));
    ifile.read((char*)&rho_Nit_avg,sizeof(double));
    ifile.read((char*)&post_proc,sizeof(bool));
    
    ifile.close();
}

void MonteCarlo::read_statistics(string filename, Array& avg, Array& std)
{
    unsigned int N = param.problem_size;
    avg.resize(N);
    std.resize(N);
    filename = filename+string("_statistics.bin");
    ifstream ifile(filename, ios::in | ios::binary);
    
    for(unsigned int i=0;i<N;i++)
    {
        ifile.read((char*)&avg(i),sizeof(double));
        ifile.read((char*)&std(i),sizeof(double));
    }
    
    ifile.close();
}

Real MonteCarlo::get_cpu_time()
{
    return (Real)clock() / CLOCKS_PER_SEC;
}

// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
static inline void loadBar(int i, int n, double start, double end)
{
    // Only update r times.
    int r=n/100;
    if(r==0)
        r=1;
    if ( i % r != 0 ) 
        return;
 
    // Calculate the ratio of complete-to-incomplete.
    int w = 60;
    float ratio = i/(float)n;
    int   c     = ratio * w;
    float time = (end-start)/60;
    float tend = time/ratio-time;
 
    cout<<"--------------------------------------------------------------"<<endl<<flush;
    
    cout<<"Complete: "<<(int)(ratio*100)<<" %"<<endl;
    cout<<setprecision(1);
    cout<<"Elapsed time: "<<time<<" minutes."<<endl;
    cout<<"Time to end : "<<tend<<" minutes."<<endl;   
    
    // Show the load bar.
    cout<<"[";
    for (int x=0; x<c; x++)
       printf("=");
 
    for (int x=c; x<w; x++)
       printf(" ");
    cout<<"]"<<endl;
    cout<<"--------------------------------------------------------------"<<endl<<flush;
    
    printf("\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J");
}

