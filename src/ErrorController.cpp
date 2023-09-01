#include <fstream>
#include "ErrorController.h"
#include "Ode.h"

ErrorController::ErrorController(Parameters* param_)
: param(param_)
{
    fixed_h = !param->dtadap;
    ode_controller = param->ode_contr;
    write_info = param->err_write_data;
    
    nu = 1.1;
    k1 = 1.0;
    k2 = 1.0;
    safe = 0.95;
    facmax = 2.;
    facmin = 0.1;
    
    errDorder=0.;
    
    reinit();
}

ErrorController::~ErrorController()
{

}

void ErrorController::reinit()
{
    errDnpu=0.;
    errDn=0.;
    errDnmu=0.;
    
    hnpu = 0.;
    hn=0.;
    hnmu=0.;
    hnmd=0.;
    
    consecutive_rej= 0;
}

void ErrorController::get_new_h(Real t, Real& h)
{
    if(write_info) 
        write_data(t);
    
    if(!fixed_h)
        hnpu= get_new_h_ODE();
    else    
        hnpu=hn;
    
    errDnmu=errDn;
    errDn=errDnpu;

    hnmd=hnmu;
    hnmu=hn;
    hn=hnpu;
    
    h=hnpu;
}

Real ErrorController::get_new_h_ODE()
{
    if(errDnpu<=nu) //accepted step
    {
        if(consecutive_rej==0)//last step accepted
            facmax=2.0;
        else //last step was rejected
        {
            facmax=1.0;
            consecutive_rej=0;
        }
    }
    else
        consecutive_rej++;        
        
    if(consecutive_rej>1) //if consecutive rejections we can estimate the 
    {                     //effective order of the error
        Real estimation_errDorder = log(errDnpu/errDn)/log(hn/hnmu);
        if(estimation_errDorder>0.1*errDorder //check if estimate makes sense
           && estimation_errDorder<=errDorder)
            fac = pow(1./errDnpu,1./estimation_errDorder);
        else
            fac = pow(1./errDnpu,1./errDorder);
    }
    else
        fac = pow(1./errDnpu,1./errDorder); // I controller
    
    if(ode_controller==PPI) //predictive proportional integral controller
    {
        fac2 = pow(1./errDnpu,k2/errDorder)*pow(errDn/errDnpu,k1/errDorder)*(hn/hnmu);
        fac =min(fac,fac2);
//        fac = fac2;
    }
    else if(ode_controller==PI)// proportional integral controller
        cout<<"Implement Proportional Integral controller"<<endl;
    
    fac = min(facmax,max(facmin,safe*fac));
    return fac*hn;
}

void ErrorController::update_hn(Real h)
{
    hn=h;
}

bool ErrorController::isAccepted()
{
    return (consecutive_rej==0);
}

Controller ErrorController::OdeController() const
{
    return ode_controller;
}

void ErrorController::set_errD_order(Real order)
{
    if(order==0)
        cout<<"ERROR: error order is zero."<<endl;
    
    errDorder=order;
}

Real& ErrorController::get_errD()
{
    return errDnpu;
}

void ErrorController::write_data(Real t)
{
    static int nout=1;
    
    if(nout==1)
    {
        ofstream outfile(sol_file+string("_errors.m"), ofstream::out);
        write_errors_file(outfile,nout,t);
        outfile.close();
    }
    else
    {
        ofstream outfile(sol_file+string("_errors.m"), ofstream::out | ofstream::app);
        write_errors_file(outfile,nout,t);
        outfile.close();
    }
    nout++;
}

void ErrorController::disp_info()
{
    cout<<"Err = "<<setw(6)<<setprecision(4);
    print_err(errDnpu);
}

void ErrorController::print_err(Real err)
{
    string red="\033[31;1m";
    string white="\033[0m";
    string yellow="\033[33;1m";
    
    if(err<=1.)
        cout<<err;
    else if(err<=nu)
        cout<<yellow<<err<<white;
    else
        cout<<red<<err<<white;
}

void ErrorController::write_errors_file(ofstream& outfile, int nout, Real t)
{
    Real phiDn = errDnpu/pow(hn,errDorder);    
    Real phiDnmu = errDn/pow(hnmu,errDorder);    
    Real phiDnmd = errDnmu/pow(hnmd,errDorder);    
    Real phiDnhatI=phiDnmu;
    Real phiDnhatPPI=phiDnmu*phiDnmu/phiDnmd;
    Real err;
    
    // write data at t_n+1
    outfile<<setprecision(16);
    outfile<<"tn("<<nout<<")="<<t<<";"<<endl;
    outfile<<"tnpu("<<nout<<")="<<t+hn<<";"<<endl;
    outfile<<"hn("<<nout<<")="<<hn<<";"<<endl;
    outfile<<"errDnpu("<<nout<<")="<<errDnpu<<";"<<endl; 
    outfile<<"errDI("<<nout<<")="<<phiDnhatI*pow(hn,errDorder)<<";"<<endl;
    outfile<<"errDPPI("<<nout<<")="<<phiDnhatPPI*pow(hn,errDorder)<<";"<<endl;
    outfile<<"phiDhatI("<<nout<<")="<<phiDnhatI<<";"<<endl;
    outfile<<"phiDhatPPI("<<nout<<")="<<phiDnhatPPI<<";"<<endl;
    outfile<<"relerrphiDI("<<nout<<")="<<abs((phiDnhatI-phiDn)/phiDn)<<";"<<endl;
    outfile<<"relerrphiDPPI("<<nout<<")="<<abs((phiDnhatPPI-phiDn)/phiDn)<<";"<<endl;
   
}