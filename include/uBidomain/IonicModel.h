#ifndef IONICMODEL_H
#define IONICMODEL_H

#include <utility>
#include <cmath>
#include <string>
#include <Eigen/Dense>
#include <iostream>

using namespace std;

typedef double Real;

typedef Eigen::MatrixXd Matrix;
typedef Eigen::RowVector2d Point;
typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXd Array;

namespace uBidomain
{
    
enum Trigger{INIT_VAL,I_APP};

class IonicModel
{
public:
    IonicModel(){};
    virtual ~IonicModel(){};
    
    virtual Vector f(const Real t, const Vector& y, bool Iapp_yes=true) const =0;        
    virtual void get_y0(Vector& y0) const =0;
    
    Real Iapp(Real t) const;
 
    Real get_tend() const;
    unsigned get_n_gating_vars() const;
    Trigger get_trigger();
    
protected:
    
    Trigger trigger;
    unsigned n_gating_vars;  
    Real tend;
    Real Imax;
    Real interval;
    Real duration;
    Real first_stim;
    Real last_stim;
};

class MyIonicModel: public IonicModel
{
public:
    MyIonicModel();
    virtual ~MyIonicModel(){};
    
    virtual Vector f(const Real t, const Vector& y, bool Iapp_yes=true) const;
    virtual Vector Iion(const Vector& y) const;
    virtual Vector dzdt(const Vector& y) const;
    
    virtual Vector Iion(const Vector& V, const Vector& z) const =0;
    virtual Vector dzdt(const Vector& V, const Vector& z) const =0;
    
    void get_y0(Vector& y0) const;
    virtual Vector get_z0(const Vector& V0) const =0;
    virtual Vector get_V0(unsigned N) const =0;
};
    
class Cubic: public MyIonicModel
{
public:
    Cubic();
    virtual ~Cubic(){};
    
    Vector Iion(const Vector& V, const Vector& z) const;
    Vector dzdt(const Vector& V, const Vector& z) const;
    
    Vector get_z0(const Vector& V0) const;
    Vector get_V0(unsigned N) const;
    
protected:
    Real eta0;
    Real Vth;
    Real Vpk;
};

class AlievPanfilov: public MyIonicModel
{
public:
    AlievPanfilov();
    virtual ~AlievPanfilov(){};
    
    Vector Iion(const Vector& V, const Vector& z) const;
    Vector dzdt(const Vector& V, const Vector& z) const;
    
    Vector get_z0(const Vector& V0) const;
    Vector get_V0(unsigned N) const;
    
protected:
    Real ga;
    Real a;
};

class MitchellSchaeffer: public MyIonicModel
{
public:
    MitchellSchaeffer();
    virtual ~MitchellSchaeffer(){};
    
    Vector Iion(const Vector& V, const Vector& z) const;
    Vector dzdt(const Vector& V, const Vector& z) const;
    
    Vector get_z0(const Vector& V0) const;
    Vector get_V0(unsigned N) const;
    
protected:
    Real tau_in;
    Real tau_out;
    Real tau_open;
    Real tau_close;
    Real V_gate;
};

class HodgkinHuxley: public MyIonicModel
{
public:
    HodgkinHuxley();
    virtual ~HodgkinHuxley(){};
    
    Vector Iion(const Vector& V, const Vector& z) const;
    Vector dzdt(const Vector& V, const Vector& z) const;
    
    Vector get_z0(const Vector& V0) const;
    Vector get_V0(unsigned N) const;
    
protected:
    Array am(const Array& V) const;
    Array an(const Array& V) const;
    Array ah(const Array& V) const;
    Array bm(const Array& V) const;
    Array bn(const Array& V) const;
    Array bh(const Array& V) const;
    
    Real gNa;
    Real gK;
    Real gL;
    Real vNa;
    Real vK;
    Real vL;
};

}

#endif /* IONICMODEL_H */

