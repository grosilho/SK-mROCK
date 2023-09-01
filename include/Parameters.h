#ifndef INIT_H
#define	INIT_H

#include "MainHeader.h"

class TimeIntegrator;
class OdeRungeKuttaIntegrator;
class DSdeRungeKuttaIntegrator;
class Ode;
class DSde;

class GetPot;

class Parameters
{
public:
    Parameters();
    ~Parameters();
    
    bool read_command_line(int argc, char** argv);
    
    bool initDSde(DSde*& sde);
    bool initDSdeTimeIntegrator(DSdeRungeKuttaIntegrator*& rk, DSde* sde);
    bool initDSdeIntegration(DSdeRungeKuttaIntegrator*& integr, DSde*& sde);
    
    void print_info();
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
    
    string rk_name;
    Real dt;
    unsigned int rho_freq;
    
    bool conv_test;
    unsigned int max_pow;
    unsigned int min_pow;
    
    unsigned int MCiter;
    bool continuous;
    int seed;
    bool post_proc;
    unsigned int n_bins;
    bool process_only;
    
    int problem_size;
};


#endif	/* INIT_H */

