#include "ClassicalOdeRungeKuttaIntegrators.h"
//#include "MatrixReplacement.h"

ExplicitEuler::ExplicitEuler(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_, ode_)
{
    Astable = false;
    s = 1;
    err_order = 1;
    
    reinit_integrator();
}

ExplicitEuler::~ExplicitEuler()
{
}

void ExplicitEuler::step(const Real t, const Real& h)
{
    Vector*& k1= integr[0];
        
    ode->f(t,*yn,*k1);
    *ynpu = *yn + h*(*k1);

    n_f_eval ++;
}

void ExplicitEuler::update_n_stages_and_h(Real& h)
{
    if(h>0.9*2.0/(this->eigmax))
    {
        h = 0.9*2.0/(this->eigmax);
        this->last=false;
    }
    
    s_max = 1;
    s_avg += s;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ImplicitEuler::ImplicitEuler(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_, ode_)
{
    Astable = true;
    s = 1;
    err_order = 1;
    
    Newton_tol = 1e-5;
    Newton_max_iter = 1e3;
    
    if(ode->has_dense_Jacobian())
    {
        J.resize(ode->get_system_size(),ode->get_system_size());
        I = Eigen::MatrixXd::Identity(ode->get_system_size(),ode->get_system_size());
    }
    else//use sparse matrices
    {
        spJ.resize(ode->get_system_size(),ode->get_system_size());
        spI.resize(ode->get_system_size(),ode->get_system_size());
        spI.setIdentity();
    }
    
    reinit_integrator();
}

ImplicitEuler::~ImplicitEuler()
{
}

void ImplicitEuler::step(const Real t, const Real& h)
{
    Vector*& dk= integr[0];
    Vector*& hfyk= integr[1];
    Vector*& ddk= integr[2];
    Vector*& tmp1= integr[3];
    Vector*& tmp2= integr[4];
    Eigen::BiCGSTAB<SpMatrix, Eigen::IncompleteLUT<Real>> solver; 
//    Eigen::ColPivHouseholderQR<Matrix> dir_solver(ode->get_system_size(),ode->get_system_size());
    Eigen::PartialPivLU<Matrix> dir_solver(ode->get_system_size());
    
    *dk *=0.;
    
    bool modified_Newton = true;
    bool ok=false;
    Newton_iter=0;
    lin_solv_iter=0;
    Real etak=1.;
    Real thetak=0.;
    Real gamma=1e-2;
    Real errk,errkm1;
    do
    {
        *tmp1 = *yn+*dk;
        ode->f(t+h,*tmp1,*hfyk);
        *hfyk *= h;
        *hfyk -= *dk;
        
        if(ode->has_dense_Jacobian())
        {
            if(Newton_iter==0 || !modified_Newton)
            {
                ode->df(t+h,*tmp1,J);
                J = I - J*h;
                dir_solver.compute(J);
            }
            *ddk = dir_solver.solve(*hfyk);
            
        }
        else
        {
            if(Newton_iter==0 || !modified_Newton)
            {
                ode->df(t+h,*tmp1,spJ);              
                spJ = spI - spJ*h;
                solver.compute(spJ);
                
                if(solver.info()!=Eigen::ComputationInfo::Success)
                {
                    cout<<"Error in iterative solver"<<endl;
                    return;
                }
            }
            *ddk = solver.solve(*hfyk);
//            *ddk = solver.solveWithGuess(*hfyk,Eigen::VectorXd(INIT TO ZERO));
            
            if(solver.info()!=Eigen::ComputationInfo::Success)
            {
                cout<<"Error in solving linear system"<<endl;
                return;
            }
            lin_solv_iter += solver.iterations();
            //cout<<"System solved with "<<solver.iterations()<<" iterations and "<<solver.error()<<" estimated error"<<endl;
        }
        
        errkm1=errk;
        errk = ddk->norm();
        
        if(Newton_iter>0)
        {
            thetak = errk/errkm1;
            etak = thetak/(1.-thetak);
        }

        Newton_iter++;
        
        if(thetak>1.)// && Newton_iter>=3)// here we detect oscillations in Newton
            break;
        
        if(etak*errk<gamma*dk->norm()*Newton_tol)
//        if(etak*errk<gamma*dk->norm()*h)
            ok=true;
        
        *dk += *ddk;

    }while(!ok && Newton_iter<Newton_max_iter);

    if(Newton_iter==Newton_max_iter && !ok)
        cout<<"WARNING: Newton algorithm did not converge: max number of iter reached."<<endl;
    else if(!ok)
        cout<<"WARNING: Newton algorithm did not converge: oscillations detected."<<endl;
    
    *ynpu = *yn+*dk;
    
    n_f_eval+=Newton_iter;
}

void ImplicitEuler::update_n_stages_and_h(Real& h)
{
    s_max = 1;
    s_avg += s;
}

void ImplicitEuler::disp_step_info(Real& t, Real& h, bool accepted)
{
    string rho =u8"\u03C1";
    cout<<setprecision(4)<<scientific;
    
    int stages = s+((param->rk_name=="ROCK2"||param->rk_name=="DROCK2"||param->rk_name=="SROCK2") ? 2:0);
    std::string delta = u8"\u0394";

    cout << scientific;
    
    cout<<"Step t = "<<setw(6)<<setprecision(4)<<t<<", "<<delta<<"t = "<<setw(8)<<setprecision(6)<<h
    <<", s = "<<setw(3)<<stages<<", "<<"Newton iter = "<<setw(4)<<Newton_iter<<", "<<"Lin solv iter = "<<setw(4)<<lin_solv_iter
    <<" and |y_n+1| = "<<setw(7)<<setprecision(4)<<ynpu->lpNorm<Eigen::Infinity>()<<". ";
    if(accepted)
        cout<<" Accepted "; 
    else
        cout<<" Rejected ";
    err_control.disp_info();
    cout<<endl;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExplicitMidpoint::ExplicitMidpoint(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_, ode_)
{
    Astable = false;
    s = 2;
    err_order = 2;
    
    reinit_integrator();
}

ExplicitMidpoint::~ExplicitMidpoint()
{
}

void ExplicitMidpoint::step(const Real t, const Real& h)
{
    Vector*& k1= integr[0];
    Vector*& k2= integr[1];
        
    ode->f(t,*yn,*k1);
    *k1 = *yn + h*(*k1)/2.;
    ode->f(t,*k1,*k2);
    *ynpu = *yn + h*(*k2);
    
    n_f_eval+=2;
}

void ExplicitMidpoint::update_n_stages_and_h(Real& h)
{
    if(h>0.9*2.0/(this->eigmax))
    {
        h = 0.9*2.0/(this->eigmax);
        this->last=false;
    }
    
    s_max = 2;
    s_avg += s;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ImplicitMidpoint::ImplicitMidpoint(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_, ode_)
{
    Astable = true;
    s = 1;
    err_order = 2;
    
    Newton_tol = 1e-8;
    Newton_max_iter = 1e3;
    
    if(ode->has_dense_Jacobian())
    {
        J.resize(ode->get_system_size(),ode->get_system_size());
        I = Eigen::MatrixXd::Identity(ode->get_system_size(),ode->get_system_size());
    }
    else//use sparse matrices
    {
        spJ.resize(ode->get_system_size(),ode->get_system_size());
        spI.resize(ode->get_system_size(),ode->get_system_size());
        spI.setIdentity();
    }
    
    reinit_integrator();
}

ImplicitMidpoint::~ImplicitMidpoint()
{
}

void ImplicitMidpoint::step(const Real t, const Real& h)
{
    Vector*& zk= integr[0];
    Vector*& hfzk= integr[1];
    Vector*& dzk= integr[2];
    
    *zk *= 0.;
    
    Eigen::BiCGSTAB<SpMatrix, Eigen::IncompleteLUT<Real>> solver; 
    Eigen::ColPivHouseholderQR<Matrix> dir_solver(ode->get_system_size(),ode->get_system_size());
    
    bool modified_Newton = true;
    bool ok=false;
    unsigned int Newton_iter=0;
    lin_solv_iter=0;
    do
    {
        *dzk = *yn + 0.5*(*zk);
        ode->f(t+h/2.,*dzk,*hfzk);
        *hfzk *= h;
        *hfzk -= *zk;
        
        if(ode->has_dense_Jacobian())
        {
            if(Newton_iter==0 || !modified_Newton)
            {
                ode->df(t+h/2.,*dzk,J);
                J = I - J*h/2.;
                dir_solver.compute(J);
            }

            *dzk = dir_solver.solve(*hfzk); 
        }
        else
        {
            if(Newton_iter==0 || !modified_Newton)
            {
                ode->df(t+h/2.,*dzk,spJ);
                spJ = spI - spJ*h/2.;
                solver.compute(spJ);
                
                if(solver.info()!=Eigen::ComputationInfo::Success)
                {
                    cout<<"Error in iterative solver"<<endl;
                    return;
                }
            }
            
            *dzk = solver.solveWithGuess(*hfzk,*zk);
            
            if(solver.info()!=Eigen::ComputationInfo::Success)
            {
                cout<<"Error in solving linear system"<<endl;
                return;
            }
            lin_solv_iter += solver.iterations();
        }
        

        
        if(dzk->norm()<zk->norm()*Newton_tol)
            ok=true;
        
        *zk += *dzk;

        Newton_iter++;

    }while(!ok && Newton_iter<Newton_max_iter);

    if(Newton_iter==Newton_max_iter && !ok)
        cout<<"WARNING: Newton algorithm did not converge."<<endl;

    *ynpu = (*yn) + (*zk);
    
    n_f_eval+=Newton_iter;
}

void ImplicitMidpoint::update_n_stages_and_h(Real& h)
{    
    s_max = 1;
    s_avg += s;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RungeKutta4::RungeKutta4(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_, ode_)
{
    Astable = false;
    s = 4;
    err_order = 4;
    
    reinit_integrator();
}

RungeKutta4::~RungeKutta4()
{
}

void RungeKutta4::step(const Real t, const Real& h)
{
    
    Vector*& k1= integr[0];
    Vector*& k2= integr[1];
    Vector*& k3= integr[2];
    Vector*& k4= integr[3];
    Vector*& tmp= integr[4];    
        
    ode->f(t,*yn,*k1);
    
    *tmp = *yn;
    *tmp+= h/2.*(*k1);
    ode->f(t+h/2.,*tmp,*k2);
    
    *tmp = *yn;
    *tmp+= h/2.*(*k2);
    ode->f(t+h/2.,*tmp,*k3);
    
    *tmp = *yn;
    *tmp+= h*(*k3);
    ode->f(t+h,*tmp,*k4);
    
    *ynpu = *yn;
    *ynpu += h/6.*(*k1);
    *ynpu += h/3.*(*k2);
    *ynpu += h/3.*(*k3);
    *ynpu += h/6.*(*k4);
    
    n_f_eval += 4;
}

void RungeKutta4::update_n_stages_and_h(Real& h)
{
    if(h>0.9*2.75/(this->eigmax))
    {
        h = 0.9*2.75/(this->eigmax);
        this->last=false;
    }
    
    s_max = 4;
    s_avg += s;
}


