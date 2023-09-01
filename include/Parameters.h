#ifndef INIT_H
#define	INIT_H

#include "MainHeader.h"

class TimeIntegrator;
class OdeRungeKuttaIntegrator;
class DSdeRungeKuttaIntegrator;
class Ode;
class DSde;

//#include "FastSlowOdeProblems.h"
//#include "SdeProblems.h"
//#include "OdeStabilizedIntegrators.h"
//#include "SdeStabilizedIntegrators.h"
//#include "TraditionalSdeRungeKuttaIntegrators.h"
//#include "ErrorController.h"

class GetPot;

class Parameters
{
public:
    Parameters();
    ~Parameters();
    
    bool read_command_line(int argc, char** argv);
    
    bool initOde(Ode*& ode);
    bool initOdeTimeIntegrator(OdeRungeKuttaIntegrator*& rk, Ode* ode);
    bool initOdeIntegration(OdeRungeKuttaIntegrator*& rk, Ode*& ode);
    
    bool initDSde(DSde*& sde);
    bool initDSdeTimeIntegrator(DSdeRungeKuttaIntegrator*& rk, DSde* sde);
    bool initDSdeIntegration(DSdeRungeKuttaIntegrator*& integr, DSde*& sde);
    
    void print_info();
    void print_info(Ode* ode);
    void print_info(DSde* sde);
    
public:
    int ntest;
    string output_path;
    string refsol_path;
    string output_file;
    string refsol_file;
    int output_freq;
    bool verbose;
    
    bool matlab_output;
    bool bin_output;
    bool specific_output;
    
    Equation eq;
    string rk_name;
    Real dt;
    unsigned int rho_freq;
    
    bool conv_test;
    unsigned int max_pow;
    unsigned int min_pow;
    
    bool dtadap;
    Real rtol;
    Real atol;
    Controller ode_contr;
    bool err_write_data;
    
    unsigned int MCiter;
    bool continuous;
    int seed;
    bool post_proc;
    unsigned int n_bins;
    bool process_only;
    
    bool parareal;
    unsigned int n_threads;
    string outer_rk_name;
    string inner_rk_name;
    Real outer_dt;
    Real inner_dt;
    
    int problem_size;
    
    //Problem specific parameters
    //Problem 20
    Real P20_dG;
    Real P20_cell_l;
    Real P20_cell_w;
    unsigned P20_nx;
    unsigned P20_ny;
    Real P20_Rl;
    Real P20_Rt;
    Real P20_si;
    Real P20_se;
    Real P20_tend;
    Real P20_Imax;
    Real P20_vertamp;
    unsigned P20_vertfreq;
    bool P20_vsmoothwave;
    Real P20_horamp;
    Real P20_horlen;
    Real P20_horpos;
    int P20_horalt;
    bool P20_hsmoothwave;
    bool P20_onlyperp;
    unsigned P20_nic;
    Real P20_horprob;
};


#endif	/* INIT_H */

