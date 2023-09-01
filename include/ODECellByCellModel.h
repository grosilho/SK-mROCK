#ifndef ODECELLBYCELLMODEL_H
#define ODECELLBYCELLMODEL_H

#include "Ode.h"
#include "CellyByCellModel.h"
#include "MyoKitIonicModel.h"
#include "Parameters.h"

class ODECellByCellModel: public virtual Ode
{
public:
    ODECellByCellModel(const int n, Parameters& param_);
    virtual ~ODECellByCellModel();
    
    void set_initial_value(Vector& y0);
    void f(Real t, Vector& x, Vector& fx);    
    
    void write_solution(const int nout, const string solname, 
                        const Real t, const Vector& y);
    
protected:
    uBidomain::CellByCellModel cbc;
    unique_ptr<uBidomain::IonicModel> im;
    Vector stim_vec;
    
    unsigned n_transmembrane_dofs;
    unsigned n_gating;
    Real Cm;
};


class MultirateODECellByCellModel: public ODECellByCellModel,
                                   public MultirateOde
{
public:
    MultirateODECellByCellModel(const int n, Parameters& param_);
    virtual ~MultirateODECellByCellModel();
    
    void fF(Real t, Vector& y, Vector& fy);
    void fS(Real t, Vector& y, Vector& fy);
    void dfF(Real t, Vector& x, Matrix& dfx);
};


#endif /* ODECELLBYCELLMODEL_H */

