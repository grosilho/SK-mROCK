#ifndef ODE_H
#define	ODE_H

#include "MainHeader.h"


class Ode
{
public:
    Ode();
    virtual ~Ode();
        
    virtual void set_initial_value(Vector& y0) = 0;    
    virtual void f(Real t, Vector& x, Vector& fx) =0;
    Real get_tend();
    
    void df(Real t, Vector& x, Matrix& dfx);  
    void df(Real t, Vector& x, SpMatrix& dfx);
    virtual void rho(Real t, Vector& y, Real& eigmax);
    
    int get_system_size();
    bool is_rho_constant();
    bool estimation_rho_known();
    bool analytical_df_known();
    bool autodiff_implemented();
    string get_problem_name();
    bool is_multirate();
    bool has_dense_Jacobian();
    virtual void write_solution(const int nout, const string solname, 
                                const Real t, const Vector& y);
    
protected:
    void FD_df(Real t, Vector& x, Matrix& fx);
    void FD_df(Real t, Vector& x, SpMatrix& fx);
    virtual void AN_df(Real t, Vector& x, Matrix& dfx);
    virtual void AD_df(Real t, Vector& x, Matrix& fx);
    virtual void AN_df(Real t, Vector& x, SpMatrix& dfx);
    virtual void AD_df(Real t, Vector& x, SpMatrix& fx);
    
    Real tend;
    
    static int neqn;
    bool cte_rho;
    bool know_rho;
    bool analytical_df;
    bool autodiff_df;
    
    string problem_name;
    
    bool multirate;
    bool dense_Jacobian;
};


// MULTIRATE ODE ---------------------------------------------------------------

class MultirateOde: public virtual Ode
{
public:
    MultirateOde();
    virtual ~MultirateOde();
    
    virtual void fF(Real t, Vector& x, Vector& fx) =0;
    virtual void fS(Real t, Vector& x, Vector& fx) =0;
    
    virtual void dfF(Real t, Vector& x, Matrix& dfx);  
    virtual void dfF(Real t, Vector& x, SpMatrix& dfx);
    
    virtual void rho(Real t, Vector& y, Real& eigmax_F, Real& eigmax_S);
};

#endif	/* ODE_H */

