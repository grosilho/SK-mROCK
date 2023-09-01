#include "ClassicalOdeRungeKuttaIntegrators.h"
//#include "MatrixReplacement.h"

ExplicitEuler::ExplicitEuler(Parameters* param_, Ode* ode_)
:OdeRungeKuttaIntegrator(param_, ode_)
{
    Astable = false;
    s = 1;
    
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
    cout<<endl;
}
