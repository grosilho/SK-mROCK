#ifndef ERRORCONTROLLER_H
#define	ERRORCONTROLLER_H

#include "MainHeader.h"
#include "Parameters.h"

class ErrorController
{
public:
    ErrorController(Parameters* param_);
    virtual ~ErrorController();
   
    void reinit();
   
    void get_new_h(Real t, Real& h);
    bool isAccepted();
    Controller OdeController() const;
    
    void update_hn(Real h);
    void set_errD_order(Real order);
    void setPhinpuhatF(Real phinpuhat);
    void setExactEDnpu(Real err, Real h);
    Real& get_errD();
    
    void disp_info();
    
protected:
    Real get_new_h_ODE();
    void print_err(Real err);
    void write_data(Real t);
    void write_errors_file(ofstream& outfile, int nout, Real t);
    
protected:
    Parameters* param;
    
    Controller ode_controller;
    bool fixed_h;
        
    Real errDnpu;
    Real errDn;
    Real errDnmu;

    Real errDorder;
    Real k1;
    Real k2;
    
    Real hnpu;
    Real hn;
    Real hnmu;
    Real hnmd;
    
    Real nu;
    Real facmax;
    Real facmin;
    Real safe;
    Real fac1;
    Real fac2;
    Real fac;
    
    int consecutive_rej;
    
    bool write_info;
    string sol_file;
};

#endif	/* ERRORCONTROLLER_H */

