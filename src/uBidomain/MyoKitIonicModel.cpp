#include "MyoKitIonicModel.h"

namespace uBidomain{    
namespace MyoKit{
            
    #define N_Vector Vector&
    #define NV_Ith_S(y, i) y.segment(i*y.size()/N_STATE,y.size()/N_STATE).array()
    using realtype = Real;
    
    #include "myokit/HodgkinHuxley1952.c" // 4  state vars
//    #include "myokit/LuoRudy1991.c"       // 8  state vars
//    #include "myokit/Fox2002.c"       // 13  state vars
//    #include "myokit/Courtemanche1998.c"  // 21 state vars, human atrial 
//    #include "myokit/OHaraRudy2011.c"     // 41 state vars

MyoKitIonicModel::MyoKitIonicModel()
{
    init();
}

MyoKitIonicModel::MyoKitIonicModel(const Parameters& param_)
{
    init();    
    
    if(param_.P20_Imax>0.)
        Imax = param_.P20_Imax;
    
    if(param_.P20_tend>0.)
        tend = param_.P20_tend;
}

void MyoKitIonicModel::init()
{
    n_gating_vars = N_STATE-1;
    
//    trigger = INIT_VAL;
    trigger = I_APP; 
    
    // For HodgkinHuxley1952
//    tend = 20.;
//    first_stim = 0.;
//    last_stim = 200.;
//    interval = 20.;
//    duration = 1.;
//    Imax = 20.; 
    
    // For LuoRudy1991, Courtemanche1998, OHaraRudy2011
    tend = 10.;
    first_stim = 0.;
    last_stim = 2000.;
    interval = 200.;
    duration = 1.;
    Imax = 120.;  //80 for cell, 250 for tissue (or more, depends)
    
    //For transversal 20x1 cells 200 seems ok
    //For longitudinal 1x20 cells  also 200
    
    updateConstants();
}

void MyoKitIonicModel::get_y0(Vector& y0) const
{
    default_initial_values(y0);
}

Vector MyoKitIonicModel::f(const Real t, const Vector& y, bool Iapp_yes) const
{
    Vector ydot(y.size());
    const unsigned N = y.size()/(1+n_gating_vars);
    
    rhs(t, y, ydot, nullptr);
    
    if(Iapp_yes)
        ydot.head(N) += Iapp(t)*Vector::Ones(N);
    
    return ydot;
}

}
}