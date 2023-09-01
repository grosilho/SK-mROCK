#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "Parameters.h"
#include "DSde.h"
#include "DSdeRungeKuttaIntegrator.h"

class MonteCarlo
{
public:
    MonteCarlo(Parameters param_);
    ~MonteCarlo();
    
    void integrate();
    void convergence_test();
    
    void print_statistics();
    
    void import_process_data(); //here no integrations are performed. Two sets of data are readen from file and processed/compared.
   
protected:
    void process_data(vector<bool>& succeed, vector<Real>& dt, 
                      vector<Real>& n_f_eval, vector<Real>& rho_Nit);
    void print_data();
    void write_data(vector<bool>& succeed);
    void compute_print_errors();
    
    void import_data(string filename, vector<Vector>& Y);
    void read_solver_info(string filename, Real& time, unsigned int& iter,
                 unsigned int& succeed_iter, string& solver_name_str,
                 Real& dt_avg, Real& n_f_eval_avg, Real& rho_Nit_avg, bool& post_proc);
    void read_statistics(string filename, Array& avg, Array& std);
    
    void print_solver_info(Real time, unsigned int iter, unsigned int succ_iter, string solver_name,
                           Real dt_avg, Real n_f_eval_avg, Real rho_Nit_avg, bool postproc);
    void print_statistics(const Array& avg, const Array& std);
    
    void compute_density_distance_area(vector<unsigned int>& n_bins, Array& dda,
                                       vector<Vector>& X1, vector<Vector>& X2);
    void estimate_self_distances(vector<unsigned int>& n_bins, unsigned int N1, unsigned int N2,
                                 Array& sd1, Array& sd2);
    
    Real get_cpu_time();
    
    Parameters param;
    vector<Vector> X;
    Real dt_avg;
    Real n_f_eval_avg;
    Real rho_Nit_avg;
    
    Array mean;
    Array std_dev;
    unsigned int succeed_MCiter;
    
    Real elapsed_time;
};

#endif /* MONTECARLO_H */
