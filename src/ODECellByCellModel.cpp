#include "ODECellByCellModel.h"

ODECellByCellModel::ODECellByCellModel(const int n, Parameters& param_)
:cbc(n,param_)
{
    cout<<"In ODECellByCellModel()"<<endl;
    
//    im = make_unique<uBidomain::AlievPanfilov>();
    im = make_unique<uBidomain::MyoKit::MyoKitIonicModel>(param_);
    
    problem_name = "ODECellByCellModel";
    
//    tend = 20.; // ms
    tend = im->get_tend();
    
    cte_rho = false;
    know_rho = false;
    analytical_df=false;
    dense_Jacobian=true;
    
    n_transmembrane_dofs = cbc.get_n_transmembrane_dofs();
    n_gating = im->get_n_gating_vars();
    neqn = (1+n_gating)*n_transmembrane_dofs;
    
    Cm = cbc.get_Cm();
    
    vector<unsigned> stim_doms={1};
    if(param_.P20_nic==2)
        stim_doms.push_back((cbc.get_N_domains()-1)/2+1);
    
    Vector global_stim_vec = Vector::Zero(cbc.get_n_dofs());
    for(auto stim_dom:stim_doms)
    {
        unsigned stim_dom_size = cbc.get_domain(stim_dom).get_tot_n_dofs();
        global_stim_vec += cbc.local_to_global_left(Vector::Ones(stim_dom_size),stim_dom,false);
    }
    stim_vec = cbc.global_to_local_left(global_stim_vec,0,false);
    
    cout<<"End of ODECellByCellModel()"<<endl;
}

ODECellByCellModel::~ODECellByCellModel()
{
}

void ODECellByCellModel::set_initial_value(Vector& y0)
{
    y0.resize(neqn);
    im->get_y0(y0);      
        
    if(im->get_trigger()==uBidomain::Trigger::INIT_VAL)
        y0.head(n_transmembrane_dofs) 
                += -stim_vec*y0(0)+stim_vec*0.5;
}

void ODECellByCellModel::f(Real t, Vector& x, Vector& fx)
{
    const auto& Vn = x.head(n_transmembrane_dofs);
    const auto& zn = x.tail(n_gating*n_transmembrane_dofs);
    
    fx = im->f(t,x,false);
    fx.head(n_transmembrane_dofs) += cbc.psi(Vn)+im->Iapp(t)*stim_vec;
    
    if(Cm!=1.0)
        fx.head(n_transmembrane_dofs) /= Cm;
}

void ODECellByCellModel::write_solution(const int nout, const string solname, 
                        const Real t, const Vector& y)
{
    if(nout==1)
    {
        ofstream out(solname+string("_geo.m"), ios::out);
        cbc.write(out,"mdom");
        out.close();
        
        cbc.mesh();
    }
    
   cbc.write_sol_and_mesh(solname, nout, y.head(n_transmembrane_dofs));
    
    
    // similar to TimeIntegrator::output_solution but without gating variables
    ofstream outfile;
    if(nout==1)
        outfile.open(solname+string("_V_evolution.bin"), ios::out | ios::binary);
    else
        outfile.open(solname+string("_V_evolution.bin"), ios::out | ios::binary | ios::app);
    outfile.write((char*)&t, sizeof(double));
    outfile.write((char*)&(y(0)), n_transmembrane_dofs*sizeof(double));
    outfile.close();
}

// -------------------------------------------------------
MultirateODECellByCellModel::MultirateODECellByCellModel(const int n, Parameters& param_)
:Ode(),
 ODECellByCellModel(n,param_),
 MultirateOde()
{
    know_rho = false;
}

MultirateODECellByCellModel::~MultirateODECellByCellModel()
{}

void MultirateODECellByCellModel::fF(Real t, Vector& x, Vector& fx)
{    
    const auto& Vn = x.head(n_transmembrane_dofs);
    
    fx.head(n_transmembrane_dofs) = cbc.psi(Vn)/Cm;
    fx.tail(n_gating*n_transmembrane_dofs) *= 0.;
}

void MultirateODECellByCellModel::fS(Real t, Vector& x, Vector& fx)
{   
    fx = im->f(t,x,false);    
    fx.head(n_transmembrane_dofs) += im->Iapp(t)*stim_vec;
    
    if(Cm!=1.0)
        fx.head(n_transmembrane_dofs) /= Cm;
}

void MultirateODECellByCellModel::dfF(Real t, Vector& x, Matrix& dfx)
{
    dfx = Matrix::Zero(n_transmembrane_dofs*(n_gating+1),n_transmembrane_dofs*(n_gating+1));
    dfx.block(0,0,n_transmembrane_dofs,n_transmembrane_dofs) = cbc.get_VtoLambdaM_nomeasure()/Cm;
}