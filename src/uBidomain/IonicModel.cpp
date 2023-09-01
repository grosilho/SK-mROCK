#include "IonicModel.h"

namespace uBidomain{

Real IonicModel::Iapp(Real t) const
{    
    if(trigger!=I_APP)
        return 0.;
       
    Real t_off = first_stim+floor((t-first_stim)/interval)*interval;
    
    if(t<=last_stim && t>=first_stim && t<=t_off+duration)
        return Imax;//*(t-t_off)/duration*exp(-(t-t_off-duration)/duration);
    else
        return 0.;
}

Real IonicModel::get_tend() const
{
    return tend;
}

unsigned IonicModel::get_n_gating_vars() const
{
    return n_gating_vars;
}

Trigger IonicModel::get_trigger()
{
    return trigger;
}

// -----------------------------------------------------------

MyIonicModel::MyIonicModel()
{
}

Vector MyIonicModel::f(const Real t, const Vector& y, bool Iapp_yes) const
{
    unsigned N = y.size()/(1+n_gating_vars);
    Vector rhs(y.size());
    
    rhs.head(N) = -Iion(y);
    if(Iapp_yes)
        rhs.head(N) += Iapp(t)*Vector::Ones(N);
    rhs.tail(N*n_gating_vars) = dzdt(y);
    
    return rhs;
}

Vector MyIonicModel::Iion(const Vector& y) const
{
    unsigned N = y.size()/(1+n_gating_vars);
    const auto& V = y.head(N);
    const auto& z = y.tail(n_gating_vars*N);
    
    return Iion(V,z);
}

Vector MyIonicModel::dzdt(const Vector& y) const
{
    unsigned N = y.size()/(1+n_gating_vars);
    const auto& V = y.head(N);
    const auto& z = y.tail(n_gating_vars*N);
    
    return dzdt(V,z);
}

void MyIonicModel::get_y0(Vector& y0) const
{
    unsigned N = y0.size()/(1+n_gating_vars);
    y0.head(N) = get_V0(N);
    y0.tail(N*n_gating_vars) = get_z0(y0.head(N));
}

// -----------------------------------------------------------

Cubic::Cubic()
{
    n_gating_vars = 0;
    
    eta0 = 10.;
    Vth = 13.;
    Vpk = 80.;
    
    trigger = INIT_VAL;
    
    tend = 5.;
    first_stim = 20.;
    last_stim = 100.;
    interval = 20.;
    duration = 1.;
    Imax = 10.; 
}

Vector Cubic::Iion(const Vector& V, const Vector& z) const
{
    return eta0*V.array()*(1.-V.array()/Vth)*(1.-V.array()/Vpk);   
}

Vector Cubic::dzdt(const Vector& V, const Vector& z) const
{
    return Vector();
}

Vector Cubic::get_z0(const Vector& V0) const
{
    return Vector();
}

Vector Cubic::get_V0(unsigned N) const
{
    return 81.*Vector::Ones(N);
}

// -----------------------------------------------------------

AlievPanfilov::AlievPanfilov()
{
    n_gating_vars = 0;
    
    ga = 8.;
    a = 0.1;
    
    trigger = INIT_VAL;
    
    tend = 5.;
    first_stim = 20.;
    last_stim = 100.;
    interval = 20.;
    duration = 1.;
    Imax = 10.; 
}

Vector AlievPanfilov::Iion(const Vector& V, const Vector& z) const
{
    return ga*V.array()*(V.array()-a)*(V.array()-1.);   
}

Vector AlievPanfilov::dzdt(const Vector& V, const Vector& z) const
{
    return Vector();
}

Vector AlievPanfilov::get_z0(const Vector& V0) const
{
    return Vector();
}

Vector AlievPanfilov::get_V0(unsigned N) const
{
    return 0.*Vector::Ones(N);
}

//--------------------------------------------------------------

MitchellSchaeffer::MitchellSchaeffer()
{
    n_gating_vars = 1;
    
//    trigger = I_APP;
    trigger = INIT_VAL;
    
    tau_in =  0.3;
    tau_out = 6.;
    tau_open = 120.;
    tau_close = 80.; // 80 or 150
    V_gate = 0.13;
    
    tend=100.;
    first_stim = 0.;
    last_stim = 2000.;
    interval = 250.;
    duration = 2.;
    Imax = 0.1;  
}

Vector MitchellSchaeffer::Iion(const Vector& V, const Vector& z) const
{
    return -(z.array()*V.array()*V.array()*(1.-V.array())/tau_in
              -V.array()/tau_out);    
}

Vector MitchellSchaeffer::dzdt(const Vector& V, const Vector& z) const
{
    return (V.array()<=V_gate).select(
                (1-z.array())/tau_open,
                -z.array()/tau_close );  
}

Vector MitchellSchaeffer::get_z0(const Vector& V0) const
{
    return (V0.array()<=V_gate).select(
                Vector::Ones(V0.size()),
                Vector::Zero(V0.size()));  
}

Vector MitchellSchaeffer::get_V0(unsigned N) const
{
    if(trigger==INIT_VAL)
        return 0.99*V_gate*Vector::Ones(N);
    else
        return Vector::Zero(N);
}

//--------------------------------------------------------------

HodgkinHuxley::HodgkinHuxley()
{
    n_gating_vars = 3;
    
//    trigger = INIT_VAL;
    trigger = I_APP;
    
    gNa = 120;
    gK = 36.;
    gL = 0.3;
    vNa = 40.;
    vK = -87.;
    vL = -64.387;
    
    tend = 3.;
    first_stim = 0.;
    last_stim = 2000.;
    interval = 200.;
    duration = 1.;
    Imax = 120.;
}

Vector HodgkinHuxley::Iion(const Vector& V, const Vector& z) const
{
    unsigned N = V.size();
    auto& m = z.segment(0,N);
    auto& n = z.segment(N,N);
    auto& h = z.segment(2*N,N);
    
        
    auto tmp = gNa*m.array().pow(3)*h.array()*(V.array()-vNa)
            +gK*n.array().pow(4)*(V.array()-vK)
            +gL*(V.array()-vL);
        
    return tmp;
}

Vector HodgkinHuxley::dzdt(const Vector& V, const Vector& z) const
{
    unsigned N = V.size();
    auto& m = z.segment(0,N);
    auto& n = z.segment(N,N);
    auto& h = z.segment(2*N,N);
    
    Vector HH(3*N);
    HH.segment(0,N)   = am(V)*(1.-m.array())-bm(V)*m.array();
    HH.segment(N,N)   = an(V)*(1.-n.array())-bn(V)*n.array();
    HH.segment(2*N,N) = ah(V)*(1.-h.array())-bh(V)*h.array();   
        
    return HH;
}

Vector HodgkinHuxley::get_z0(const Vector& V0) const
{
    unsigned N = V0.size();
    Vector z0(n_gating_vars*N); //m,n,h
    
    z0.segment(0,N) = 0.05*Vector::Ones(N);
    z0.segment(N,N) = 0.317*Vector::Ones(N);
    z0.segment(2*N,N) = 0.595*Vector::Ones(N);
    
    //equilibrium of gating variables
//    z0.segment(0,N) = am(V0)/(am(V0)+bm(V0));
//    z0.segment(N,N) = an(V0)/(an(V0)+bn(V0));
//    z0.segment(2*N,N) = ah(V0)/(ah(V0)+bh(V0));
    
    return z0;
}

Vector HodgkinHuxley::get_V0(unsigned N) const
{
    if(trigger==INIT_VAL)
        return -65.*Vector::Ones(N);
    else
        return -75.*Vector::Ones(N);
}

Array HodgkinHuxley::am(const Array& V) const
{   
    return -0.1*(V+50.)/(Eigen::exp(-(V+50.)/10.)-1.);
}

Array HodgkinHuxley::an(const Array& V) const
{    
    return -0.01*(V+65.)/(Eigen::exp(-(V+65.)/10.)-1.);
}

Array HodgkinHuxley::ah(const Array& V) const
{
    return 0.07*Eigen::exp(-(V+75)/20.);
}

Array HodgkinHuxley::bm(const Array& V) const
{
    return 4.0*Eigen::exp(-(V+75)/18.);
}

Array HodgkinHuxley::bn(const Array& V) const
{
    return 0.125*Eigen::exp(-(V+75)/80.);
}

Array HodgkinHuxley::bh(const Array& V) const
{
    return 1.0/(Eigen::exp(-(V+45.)/10.)+1.);
}

}