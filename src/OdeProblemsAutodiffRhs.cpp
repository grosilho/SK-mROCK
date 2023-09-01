#include "OdeProblems.h"

//------------------------------------------------------------------------------
ADVector DahlquistTestProblem::AD_f(Real t, ADVector& x)
{   
    static ADVector fx = ADVector::Zero(x.size());
    fx(0) = (lambda+xi)*x(0);
    return fx;
}

void DahlquistTestProblem::AD_df(Real t, Vector& x, Matrix& dfx)
{
    ADVector AD_x = x;
    dfx = autodiff::jacobian(AD_f, autodiff::detail::wrt(AD_x), autodiff::detail::at(t,AD_x));
}

//------------------------------------------------------------------------------
ADVector ScalarNonStiffNonLinearTest::AD_f(Real t, ADVector& x)
{
    static ADVector fx = ADVector::Zero(x.size());
    fx(0) = 0.25*x(0)+0.5*sqrt(x(0)*x(0)+1.);
    return fx;
}

void ScalarNonStiffNonLinearTest::AD_df(Real t, Vector& x, Matrix& dfx)
{    
    ADVector AD_x = x;
    dfx = autodiff::jacobian(AD_f, autodiff::detail::wrt(AD_x), autodiff::detail::at(t,AD_x));
}

//------------------------------------------------------------------------------

ADVector NeuronCable::AD_f(Real t, ADVector& x)
{   
    static ADVector fx = ADVector::Zero(x.size());
    fx(0)=nu*2.*(x(1)-x(0))*(neqn-1.)*(neqn-1.)-beta*x(0);
    fx(neqn-1)=nu*2.*(x(neqn-2)-x(neqn-1))*(neqn-1.)*(neqn-1.)-beta*x(neqn-1);
    for (int i=1;i<neqn-1;i++)
    {
        fx(i)=nu*(x(i-1)-2.*x(i)+x(i+1))*(neqn-1.)*(neqn-1.)-beta*x(i);
        if(abs(i/(neqn-1.)-0.5)<0.1)
            fx(i) += 5.*exp(1.-1e-2/(1e-2-(i/(neqn-1.)-0-5)*(i/(neqn-1.)-0.5)));
    }
    return fx;
}

void NeuronCable::AD_df(Real t, Vector& x, SpMatrix& dfx)
{
    ADVector AD_x = x;
    dfx = autodiff::jacobian(AD_f, autodiff::detail::wrt(AD_x), autodiff::detail::at(t,AD_x)).sparseView();
}

//------------------------------------------------------------------------------

ADVector Brusselator::AD_f(Real t, ADVector& x)
{   
    static ADVector fx = ADVector::Zero(x.size());
    fx(0) = (alpha-1.)*x(0)+alpha*x(0)*x(0)+(x(0)+1.)*(x(0)+1.)*x(1);
    fx(1) = -alpha*x(0)-alpha*x(0)*x(0)-(x(0)+1.)*(x(0)+1.)*x(1);
    return fx;
}

void Brusselator::AD_df(Real t, Vector& x, Matrix& dfx)
{    
    ADVector AD_x = x;
    dfx = autodiff::jacobian(AD_f, autodiff::detail::wrt(AD_x), autodiff::detail::at(t,AD_x));
}

//------------------------------------------------------------------------------

ADVector PDEBrusselator::AD_f(Real t, ADVector& x)
{   
    static ADVector fx = ADVector::Zero(x.size());
    
    Real ul = 1.;
    Real ur = 1.;
    Real vl = 3.;
    Real vr = 3.;
    int mod,j1,j2;
    
    fx(0) = A + x(0)*x(0)*(x(Nu)+vl)/2. - (B+1.)*x(0) + alpha*(ul-2.*x(0)+x(1))/Hu/Hu;
    fx(Nu-1) = A + x(Nu-1)*x(Nu-1)*(x(Nu+Nv-1)+vr)/2. - (B+1.)*x(Nu-1) + alpha*(x(Nu-2)-2.*x(Nu-1)+ur)/Hu/Hu;
    for(int i=1;i<Nu-1;i++)
    {
        mod = i % 2;
        j2 = Nu + (i-mod)/2;
        j1 = j2-1+mod;
        fx(i) = A + x(i)*x(i)*(x(j1)+x(j2))/2. - (B+1.)*x(i) + alpha*(x(i-1)-2.*x(i)+x(i+1))/Hu/Hu;
    }
    
    fx(Nu) = B*x(1) - x(1)*x(1)*x(Nu) + alpha*(vl-2.*x(Nu)+x(Nu+1))/Hv/Hv;
    fx(Nu+Nv-1) = B*x(Nu-2) - x(Nu-2)*x(Nu-2)*x(Nu+Nv-1) + alpha*(x(Nu+Nv-2)-2.*x(Nu+Nv-1)+vr)/Hv/Hv;
    for(int i=Nu+1;i<Nu+Nv-1;i++)
    {
        j1 = 2*(i-Nu)+1;
        fx(i) = B*x(j1) - x(j1)*x(j1)*x(i) + alpha*(x(i-1)-2.*x(i)+x(i+1))/Hv/Hv;
    }
    
    return fx;
}

void PDEBrusselator::AD_df(Real t, Vector& x, SpMatrix& dfx)
{    
    ADVector AD_x = x;
    dfx = autodiff::jacobian(AD_f, autodiff::detail::wrt(AD_x), autodiff::detail::at(t,AD_x)).sparseView();
}

//------------------------------------------------------------------------------


ADVector Krogh10::AD_f(Real t, ADVector& x)
{   
    static ADVector fx = ADVector::Zero(x.size());
 
    fx(0) = x(2);
    fx(1) = x(3);
    fx(2) = 2.*x(3)+x(0)-mus*(x(0)+mu)/pow(sqrt(((x(0)+mu)*(x(0)+mu)+x(1)*x(1))),3)-mu*(x(0)-mus)/pow(sqrt(((x(0)-mus)*(x(0)-mus)+x(1)*x(1))),3);
    fx(3) = -2.*x(2)+x(1)-mus*x(1)/pow(sqrt(((x(0)+mu)*(x(0)+mu)+x(1)*x(1))),3)-mu*x(1)/pow(sqrt(((x(0)-mus)*(x(0)-mus)+x(1)*x(1))),3);
    
    return fx;
}

void Krogh10::AD_df(Real t, Vector& x, Matrix& dfx)
{    
    ADVector AD_x = x;
    dfx = autodiff::jacobian(AD_f, autodiff::detail::wrt(AD_x), autodiff::detail::at(t,AD_x));
}

//void PopulationDynamics::f(Real t, Vector& x, Vector& fx)
//{   
//    fx(0)=alpha*(x(1)-1.)-lambda1*x(0)*(1-x(0));
//    fx(1)=-lambda2*x(1)*(1-x(1));
//}
//
//void VanDerPol::f(Real t, Vector& x, Vector& fx)
//{   
//    fx(0) = x(1);
//    fx(1) = ((1.-x(0)*x(0))*x(1)-x(0))/eps;
//}
//
//void OneDimHodgkinHuxley::f(Real t, Vector& z, Vector& fz)
//{           
//    Real V,n,m,h,Ionic;
//    
//    V = z(0);    
//    m = z(1);
//    h = z(2);
//    n = z(3);
//
//    Ionic = gNa*m*m*m*h*(V-VNa) + gK*n*n*n*n*(V-VK) + gL*(V-VL);
//    fz(0) = (Iapplied(t)-Ionic)/CM;    
//    fz(1) = (alpha_m(V)*(1.-m)-beta_m(V)*m)/t_m;
//    fz(2) = (alpha_h(V)*(1.-h)-beta_h(V)*h)/t_h;
//    fz(3) = (alpha_n(V)*(1.-n)-beta_n(V)*n)/t_n;
//}

//
//void DiffusionRefinedMesh::f(Real t, Vector& x, Vector& fx)
//{   
//    fx(0)=nu*2.*(x(1)-x(0))/H1/H1; // Neumann bnd cond
//    fx(neqn-1)=nu*2.*(x(neqn-2)-x(neqn-1))/H2/H2;
//    fx(N1) = nu*(H2*x(N1-1)-(H1+H2)*x(N1)+H1*x(N1+1))/(0.5*H1*H2*(H1+H2));
//    for(int i=1;i<N1;i++)
//        fx(i)=nu*(x(i-1)-2.*x(i)+x(i+1))/H1/H1;
//    for(int i=N1+1;i<N1+N2;i++)
//        fx(i)=nu*(x(i-1)-2.*x(i)+x(i+1))/H2/H2;
//}
//
//void InfectiousDiseaseTransmission::f(Real t, Vector& z, Vector& fz)
//{           
//    Real S,I,V, fV;
//    S = z(0);
//    I = z(1);
//    V = z(2);
//    
//    fV = 1.-exp(-a*V);
//    
//    fz(0) = p-m*S-r*S*fV;
//    fz(1) = r*S*fV-(m+g)*I;
//    fz(2) = o*I-e*V;
//    
//    if(abs(t-25.)<5.)
//        fz(2) += 1e5;
//}
//
ADVector RobertsonChemicalSystem::AD_f(Real t, ADVector& z)
{                
    static ADVector fz = ADVector::Zero(z.size());
    fz(0) = -k1*z(0) + k2*z(1)*z(2);
    fz(1) =  k1*z(0) - k2*z(1)*z(2)-k3*z(1)*z(1);
    fz(2) =                         k3*z(1)*z(1);
    return fz;
}

void RobertsonChemicalSystem::AD_df(Real t, Vector& x, Matrix& dfx)
{    
    ADVector AD_x = x;
    dfx = autodiff::jacobian(AD_f, autodiff::detail::wrt(AD_x), autodiff::detail::at(t,AD_x));
}
//
//void Oregonator::f(Real t, Vector& x, Vector& fx)
//{   
//    fx(0) = -k1*x(0)*x(1) - k3*x(0)*x(2);
//    fx(1) = -k1*x(0)*x(1) - k2*x(1)*x(2) + k5*x(4);
//    fx(2) =  k1*x(0)*x(1) - k2*x(1)*x(2) + k3*x(0)*x(2) - 2.*k4*x(2)*x(2);   
//    fx(3) =  k2*x(1)*x(2) + k4*x(2)*x(2);
//    fx(4) =  k3*x(0)*x(2) - k5*x(4);
//}
//
//void CUSP::f(Real t, Vector& x, Vector& fx)
//{   
//    int n=neqn/3;
//    Real D = sigma*n*n;
//    
//    fx(0) = D*(x(n-1)-2.*x(0)+x(1)) - (pow(x(0),3)+x(n)*x(0)+x(2*n))/eps ;
//    fx(n-1) = D*(x(n-2)-2.*x(n-1)+x(0)) - (pow(x(n-1),3)+x(2*n-1)*x(n-1)+x(3*n-1))/eps;
//    for(int i=1;i<n-1;i++)
//        fx(i) = D*(x(i-1)-2.*x(i)+x(i+1)) - (pow(x(i),3)+x(n+i)*x(i)+x(2*n+i))/eps;
//    
//    fx(n) = D*(x(2*n-1)-2.*x(n)+x(n+1)) + (x(2*n)+0.07*v(x(0)));
//    fx(2*n-1) = D*(x(2*n-2)-2.*x(2*n-1)+x(n)) + (x(3*n-1)+0.07*v(x(n-1)));
//    for(int i=n+1;i<2*n-1;i++)
//        fx(i) = D*(x(i-1)-2.*x(i)+x(i+1)) + (x(i+n)+0.07*v(x(i-n)));
//    
//    fx(2*n) = D*(x(3*n-1)-2.*x(2*n)+x(2*n+1)) + (1.-pow(x(n),2))*x(2*n)-x(n)-0.4*x(0)+0.035*v(x(0));
//    fx(3*n-1) = D*(x(3*n-2)-2.*x(3*n-1)+x(2*n)) + (1.-pow(x(2*n-1),2))*x(3*n-1)-x(2*n-1)-0.4*x(n-1)+0.035*v(x(n-1));
//    for(int i=2*n+1;i<3*n-1;i++)
//        fx(i) = D*(x(i-1)-2.*x(i)+x(i+1)) + (1.-pow(x(i-n),2))*x(i)-x(i-n)-0.4*x(i-2*n)+0.035*v(x(i-2*n));
//}
//
//void mReactionDiffusion2DEquations::f(Real t, Vector& y, Vector& fy)
//{   
//    unsigned int kNsq;
//    
//    // Intra-layer diffusion
//    for(unsigned int k=0; k<m; k++)
//    {
//        kNsq = k*Nsq;
//        fy(kNsq) =  D[k]*(y(kNsq+N-1)-2.*y(kNsq)+y(kNsq+1))/h/h                                     //Lower-Left corner
//                  + D[k]*(y(kNsq+N)-2.*y(kNsq)+y(kNsq+Nsq-N))/h/h;
//        fy(kNsq+N-1) =  D[k]*(y(kNsq)-2.*y(kNsq+N-1)+y(kNsq+N-2))/h/h                               //Lower-Right corner
//                      + D[k]*(y(kNsq+2*N-1)-2.*y(kNsq+N-1)+y(kNsq+Nsq-1))/h/h;
//        fy(kNsq+Nsq-N) =  D[k]*(y(kNsq+Nsq-1)-2.*y(kNsq+Nsq-N)+y(kNsq+Nsq-N+1))/h/h                 //Upper-Left corner
//                        + D[k]*(y(kNsq+Nsq-2*N)-2.*y(kNsq+Nsq-N)+y(kNsq))/h/h;
//        fy(kNsq+Nsq-1) =  D[k]*(y(kNsq+Nsq-2)-2.*y(kNsq+Nsq-1)+y(kNsq+Nsq-N))/h/h                   //Upper-Right corner
//                        + D[k]*(y(kNsq+Nsq-N-1)-2.*y(kNsq+Nsq-1)+y(kNsq+N-1))/h/h;
//        for(unsigned int j=1;j<N-1;j++)
//        {
//            fy(kNsq+j) =  D[k]*(y(kNsq+j-1)-2.*y(kNsq+j)+y(kNsq+j+1))/h/h                           // Lower boundary
//                        + D[k]*(y(kNsq+j+N)-2.*y(kNsq+j)+y(kNsq+j+Nsq-N))/h/h;
//            fy(kNsq+Nsq-N+j) =  D[k]*(y(kNsq+Nsq-N+j-1)-2.*y(kNsq+Nsq-N+j)+y(kNsq+Nsq-N+j+1))/h/h   // Upper boundary
//                              + D[k]*(y(kNsq+Nsq-N+j-N)-2.*y(kNsq+Nsq-N+j)+y(kNsq+j))/h/h;
//            fy(kNsq+j*N) =  D[k]*(y(kNsq+j*N+1)-2.*y(kNsq+j*N)+y(kNsq+j*N+N-1))/h/h                 // Left boundary
//                          + D[k]*(y(kNsq+j*N-N)-2.*y(kNsq+j*N)+y(kNsq+j*N+N))/h/h;
//            fy(kNsq+j*N+N-1) =  D[k]*(y(kNsq+j*N+N-2)-2.*y(kNsq+j*N+N-1)+y(kNsq+j*N))/h/h           // Right boundary
//                              + D[k]*(y(kNsq+j*N+2*N-1)-2.*y(kNsq+j*N+N-1)+y(kNsq+j*N-1))/h/h;
//            
//            for(unsigned int i=1;i<N-1;i++)
//                fy(kNsq+j*N+i) =  D[k]*(y(kNsq+j*N+i+1)-2.*y(kNsq+j*N+i)+y(kNsq+j*N+i-1))/h/h       // Internal points
//                                + D[k]*(y(kNsq+j*N-N+i)-2.*y(kNsq+j*N+i)+y(kNsq+j*N+N+i))/h/h;
//        }        
//   }
//    
//    // Inter-layer diffusion and reaction terms
//    for(unsigned int k=0; k<m; k++)    
//    for(unsigned int i=0;i<N;i++)
//    for(unsigned int j=0;j<N;j++)
//        fy(k*Nsq + N*i + j) += inter_diff_and_reac(y, k, i,j);
//    
//}
//
//void RadiationDiffusion::f(Real t, Vector& x, Vector& fx)
//{
//    //index i*N+j corresponds to E(x) with x=(j*h,i*h)
//    //index Nsq+i*N+j corresponds to T(x) with x=(j*h,i*h)
//    Real Eij,Eijp1,Eijm1,Eip1j,Eim1j,Eip1jp1,Eim1jp1,Eip1jm1,Eim1jm1;
//    Real Tij,Tijp1,Tijm1,Tip1j,Tim1j;
//    Real Eijp05,Eijm05,Eip05j,Eim05j;
//    Real Tijp05,Tijm05,Tip05j,Tim05j;
//    Real gradEijp05,gradEijm05,gradEip05j,gradEim05j;
//    Real radiation;
//    
//    //Diffusion and radiation in the interior of the domain
//    for(unsigned int i=1;i<N-1;i++)
//    for(unsigned int j=1;j<N-1;j++)
//    {
//        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//        Eij = x(N*i+j);
//        Eijp1 = x(N*i+j+1);
//        Eijm1 = x(N*i+j-1);
//        Eip1j = x(N*(i+1)+j);
//        Eim1j = x(N*(i-1)+j);
//        Eip1jp1 = x(N*(i+1)+j+1);
//        Eim1jp1 = x(N*(i-1)+j+1);
//        Eip1jm1 = x(N*(i+1)+j-1);
//        Eim1jm1 = x(N*(i-1)+j-1);
//        Tij = x(Nsq+N*i+j);
//        Tijp1 = x(Nsq+N*i+j+1);
//        Tijm1 = x(Nsq+N*i+j-1);
//        Tip1j = x(Nsq+N*(i+1)+j);
//        Tim1j = x(Nsq+N*(i-1)+j);
//        Eijp05 = (Eij+Eijp1)/2.;
//        Eijm05 = (Eij+Eijm1)/2.;
//        Eip05j = (Eij+Eip1j)/2.;
//        Eim05j = (Eij+Eim1j)/2.;
//        Tijp05 = (Tij+Tijp1)/2.;
//        Tijm05 = (Tij+Tijm1)/2.;
//        Tip05j = (Tij+Tip1j)/2.;
//        Tim05j = (Tij+Tim1j)/2.;
//        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//        gradEijp05 = sqrt(gradEijp05);
//        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//        gradEijm05 = sqrt(gradEijm05);
//        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//        gradEip05j = sqrt(gradEip05j);
//        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//        gradEim05j = sqrt(gradEim05j);
//        
//        // computing diffusion
//        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                                -D2(Tijm05)*(Tij  -Tijm1)
//                                +D2(Tip05j)*(Tip1j-Tij)
//                                -D2(Tim05j)*(Tij  -Tim1j));
//        
//        //adding radiation
//        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//        fx(i*N+j) += radiation;
//        fx(Nsq+i*N+j) -= radiation;
//    }
//    
//    // upper boundary
//    for(unsigned int j=1;j<N-1;j++)
//    {
//        //same code as in the interior of the domain but we use homogeneous Neumann 
//        // and all i+1 become i-1
//        unsigned int i=N-1;
//        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//        Eij = x(N*i+j);
//        Eijp1 = x(N*i+j+1);
//        Eijm1 = x(N*i+j-1);
//        Eip1j = x(N*(i-1)+j);
//        Eim1j = x(N*(i-1)+j);
//        Eip1jp1 = x(N*(i-1)+j+1);
//        Eim1jp1 = x(N*(i-1)+j+1);
//        Eip1jm1 = x(N*(i-1)+j-1);
//        Eim1jm1 = x(N*(i-1)+j-1);
//        Tij = x(Nsq+N*i+j);
//        Tijp1 = x(Nsq+N*i+j+1);
//        Tijm1 = x(Nsq+N*i+j-1);
//        Tip1j = x(Nsq+N*(i-1)+j);
//        Tim1j = x(Nsq+N*(i-1)+j);
//        Eijp05 = (Eij+Eijp1)/2.;
//        Eijm05 = (Eij+Eijm1)/2.;
//        Eip05j = (Eij+Eip1j)/2.;
//        Eim05j = (Eij+Eim1j)/2.;
//        Tijp05 = (Tij+Tijp1)/2.;
//        Tijm05 = (Tij+Tijm1)/2.;
//        Tip05j = (Tij+Tip1j)/2.;
//        Tim05j = (Tij+Tim1j)/2.;
//        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//        gradEijp05 = sqrt(gradEijp05);
//        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//        gradEijm05 = sqrt(gradEijm05);
//        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//        gradEip05j = sqrt(gradEip05j);
//        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//        gradEim05j = sqrt(gradEim05j);
//        
//        // computing diffusion
//        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                                -D2(Tijm05)*(Tij  -Tijm1)
//                                +D2(Tip05j)*(Tip1j-Tij)
//                                -D2(Tim05j)*(Tij  -Tim1j));
//        
//        //adding radiation
//        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//        fx(i*N+j) += radiation;
//        fx(Nsq+i*N+j) -= radiation;
//    }
//    
//    // lower boundary
//    for(unsigned int j=1;j<N-1;j++)
//    {
//        unsigned int i=0;
//        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//        Eij = x(N*i+j);
//        Eijp1 = x(N*i+j+1);
//        Eijm1 = x(N*i+j-1);
//        Eip1j = x(N*(i+1)+j);
//        Eim1j = x(N*(i+1)+j);
//        Eip1jp1 = x(N*(i+1)+j+1);
//        Eim1jp1 = x(N*(i+1)+j+1);
//        Eip1jm1 = x(N*(i+1)+j-1);
//        Eim1jm1 = x(N*(i+1)+j-1);
//        Tij = x(Nsq+N*i+j);
//        Tijp1 = x(Nsq+N*i+j+1);
//        Tijm1 = x(Nsq+N*i+j-1);
//        Tip1j = x(Nsq+N*(i+1)+j);
//        Tim1j = x(Nsq+N*(i+1)+j);
//        Eijp05 = (Eij+Eijp1)/2.;
//        Eijm05 = (Eij+Eijm1)/2.;
//        Eip05j = (Eij+Eip1j)/2.;
//        Eim05j = (Eij+Eim1j)/2.;
//        Tijp05 = (Tij+Tijp1)/2.;
//        Tijm05 = (Tij+Tijm1)/2.;
//        Tip05j = (Tij+Tip1j)/2.;
//        Tim05j = (Tij+Tim1j)/2.;
//        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//        gradEijp05 = sqrt(gradEijp05);
//        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//        gradEijm05 = sqrt(gradEijm05);
//        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//        gradEip05j = sqrt(gradEip05j);
//        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//        gradEim05j = sqrt(gradEim05j);
//        
//        // computing diffusion
//        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                                -D2(Tijm05)*(Tij  -Tijm1)
//                                +D2(Tip05j)*(Tip1j-Tij)
//                                -D2(Tim05j)*(Tij  -Tim1j));
//        
//        //adding radiation
//        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//        fx(i*N+j) += radiation;
//        fx(Nsq+i*N+j) -= radiation;
//    }
//    
//    //left boundary
//    for(unsigned int i=1;i<N-1;i++)
//    {
//        unsigned int j=0;
//        //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=1 for E
//        //yielding that value for Eijm1
//        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//        Eij = x(N*i+j);
//        Eijp1 = x(N*i+j+1);        
//        Eip1j = x(N*(i+1)+j);
//        Eim1j = x(N*(i-1)+j);
//        Eip1jp1 = x(N*(i+1)+j+1);
//        Eim1jp1 = x(N*(i-1)+j+1);        
//        Tij = x(Nsq+N*i+j);
//        Tijp1 = x(Nsq+N*i+j+1);
//        Tijm1 = x(Nsq+N*i+j+1);
//        Tip1j = x(Nsq+N*(i+1)+j);
//        Tim1j = x(Nsq+N*(i-1)+j);
//        Eijm1 = Eijp1+3.*h*sigma(j*h,i*h,Tij)*(4.-Eij);
//        Eip1jm1 = Eip1jp1+3.*h*sigma(j*h,(i+1)*h,Tip1j)*(4.-Eip1j);
//        Eim1jm1 = Eim1jp1+3.*h*sigma(j*h,(i-1)*h,Tim1j)*(4.-Eim1j);
//        Eijp05 = (Eij+Eijp1)/2.;
//        Eijm05 = (Eij+Eijm1)/2.;
//        Eip05j = (Eij+Eip1j)/2.;
//        Eim05j = (Eij+Eim1j)/2.;
//        Tijp05 = (Tij+Tijp1)/2.;
//        Tijm05 = (Tij+Tijm1)/2.;
//        Tip05j = (Tij+Tip1j)/2.;
//        Tim05j = (Tij+Tim1j)/2.;
//        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//        gradEijp05 = sqrt(gradEijp05);
//        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//        gradEijm05 = sqrt(gradEijm05);
//        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//        gradEip05j = sqrt(gradEip05j);
//        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//        gradEim05j = sqrt(gradEim05j);
//        
//        // computing diffusion
//        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                                -D2(Tijm05)*(Tij  -Tijm1)
//                                +D2(Tip05j)*(Tip1j-Tij)
//                                -D2(Tim05j)*(Tij  -Tim1j));
//        
//        //adding radiation
//        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//        fx(i*N+j) += radiation;
//        fx(Nsq+i*N+j) -= radiation;
//    }
//    
//    //right boundary
//    for(unsigned int i=1;i<N-1;i++)
//    {
//        //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=0 for E
//        //yielding that value for Eijp1
//        unsigned int j=N-1;
//        //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//        Eij = x(N*i+j);        
//        Eijm1 = x(N*i+j-1);
//        Eip1j = x(N*(i+1)+j);
//        Eim1j = x(N*(i-1)+j);
//        Eip1jm1 = x(N*(i+1)+j-1);
//        Eim1jm1 = x(N*(i-1)+j-1);
//        Tij = x(Nsq+N*i+j);
//        Tijp1 = x(Nsq+N*i+j-1);
//        Tijm1 = x(Nsq+N*i+j-1);
//        Tip1j = x(Nsq+N*(i+1)+j);
//        Tim1j = x(Nsq+N*(i-1)+j);
//        Eijp1 = Eijm1-3.*h*sigma(j*h,i*h,Tij)*Eij;
//        Eip1jp1 = Eip1jm1-3.*h*sigma(j*h,(i+1)*h,Tip1j)*Eip1j;
//        Eim1jp1 = Eim1jm1-3.*h*sigma(j*h,(i-1)*h,Tim1j)*Eim1j;
//        Eijp05 = (Eij+Eijp1)/2.;
//        Eijm05 = (Eij+Eijm1)/2.;
//        Eip05j = (Eij+Eip1j)/2.;
//        Eim05j = (Eij+Eim1j)/2.;
//        Tijp05 = (Tij+Tijp1)/2.;
//        Tijm05 = (Tij+Tijm1)/2.;
//        Tip05j = (Tij+Tip1j)/2.;
//        Tim05j = (Tij+Tim1j)/2.;
//        gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//        gradEijp05 = sqrt(gradEijp05);
//        gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//        gradEijm05 = sqrt(gradEijm05);
//        gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//        gradEip05j = sqrt(gradEip05j);
//        gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//        gradEim05j = sqrt(gradEim05j);
//        
//        // computing diffusion
//        fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                            -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                            +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                            -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//        fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                                -D2(Tijm05)*(Tij  -Tijm1)
//                                +D2(Tip05j)*(Tip1j-Tij)
//                                -D2(Tim05j)*(Tij  -Tim1j));
//        
//        //adding radiation
//        radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//        fx(i*N+j) += radiation;
//        fx(Nsq+i*N+j) -= radiation;
//    }
//    
//    // upper left corner
//    unsigned int i=N-1;
//    unsigned int j=0;
//    //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=1 for E
//    //yielding that value for Eijm1
//    //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//    Eij = x(N*i+j);
//    Eijp1 = x(N*i+j+1);        
//    Eip1j = x(N*(i-1)+j);
//    Eim1j = x(N*(i-1)+j);
//    Eip1jp1 = x(N*(i-1)+j+1);
//    Eim1jp1 = x(N*(i-1)+j+1);        
//    Tij = x(Nsq+N*i+j);
//    Tijp1 = x(Nsq+N*i+j+1);
//    Tijm1 = x(Nsq+N*i+j+1);
//    Tip1j = x(Nsq+N*(i-1)+j);
//    Tim1j = x(Nsq+N*(i-1)+j);
//    Eijm1 = Eijp1+3.*h*sigma(j*h,i*h,Tij)*(4.-Eij);
//    Eip1jm1 = Eip1jp1+3.*h*sigma(j*h,(i+1)*h,Tip1j)*(4.-Eip1j);
//    Eim1jm1 = Eim1jp1+3.*h*sigma(j*h,(i-1)*h,Tim1j)*(4.-Eim1j);
//    Eijp05 = (Eij+Eijp1)/2.;
//    Eijm05 = (Eij+Eijm1)/2.;
//    Eip05j = (Eij+Eip1j)/2.;
//    Eim05j = (Eij+Eim1j)/2.;
//    Tijp05 = (Tij+Tijp1)/2.;
//    Tijm05 = (Tij+Tijm1)/2.;
//    Tip05j = (Tij+Tip1j)/2.;
//    Tim05j = (Tij+Tim1j)/2.;
//    gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//    gradEijp05 = sqrt(gradEijp05);
//    gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//    gradEijm05 = sqrt(gradEijm05);
//    gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//    gradEip05j = sqrt(gradEip05j);
//    gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//    gradEim05j = sqrt(gradEim05j);
//
//    // computing diffusion
//    fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                        -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                        +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                        -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//    fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                            -D2(Tijm05)*(Tij  -Tijm1)
//                            +D2(Tip05j)*(Tip1j-Tij)
//                            -D2(Tim05j)*(Tij  -Tim1j));
//
//    //adding radiation
//    radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//    fx(i*N+j) += radiation;
//    fx(Nsq+i*N+j) -= radiation;
//    
//    
//    //lower left corner
//    i=0;
//    j=0;
//    //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=1 for E
//    //yielding that value for Eijm1
//    //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//    Eij = x(N*i+j);
//    Eijp1 = x(N*i+j+1);        
//    Eip1j = x(N*(i+1)+j);
//    Eim1j = x(N*(i+1)+j);
//    Eip1jp1 = x(N*(i+1)+j+1);
//    Eim1jp1 = x(N*(i+1)+j+1);        
//    Tij = x(Nsq+N*i+j);
//    Tijp1 = x(Nsq+N*i+j+1);
//    Tijm1 = x(Nsq+N*i+j+1);
//    Tip1j = x(Nsq+N*(i+1)+j);
//    Tim1j = x(Nsq+N*(i+1)+j);
//    Eijm1 = Eijp1+3.*h*sigma(j*h,i*h,Tij)*(4.-Eij);
//    Eip1jm1 = Eip1jp1+3.*h*sigma(j*h,(i+1)*h,Tip1j)*(4.-Eip1j);
//    Eim1jm1 = Eim1jp1+3.*h*sigma(j*h,(i-1)*h,Tim1j)*(4.-Eim1j);
//    Eijp05 = (Eij+Eijp1)/2.;
//    Eijm05 = (Eij+Eijm1)/2.;
//    Eip05j = (Eij+Eip1j)/2.;
//    Eim05j = (Eij+Eim1j)/2.;
//    Tijp05 = (Tij+Tijp1)/2.;
//    Tijm05 = (Tij+Tijm1)/2.;
//    Tip05j = (Tij+Tip1j)/2.;
//    Tim05j = (Tij+Tim1j)/2.;
//    gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//    gradEijp05 = sqrt(gradEijp05);
//    gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//    gradEijm05 = sqrt(gradEijm05);
//    gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//    gradEip05j = sqrt(gradEip05j);
//    gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//    gradEim05j = sqrt(gradEim05j);
//
//    // computing diffusion
//    fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                        -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                        +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                        -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//    fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                            -D2(Tijm05)*(Tij  -Tijm1)
//                            +D2(Tip05j)*(Tip1j-Tij)
//                            -D2(Tim05j)*(Tij  -Tim1j));
//
//    //adding radiation
//    radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//    fx(i*N+j) += radiation;
//    fx(Nsq+i*N+j) -= radiation;
//    
//    
//    //upper right corner
//    //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=0 for E
//    //yielding that value for Eijp1
//    i=N-1;
//    j=N-1;
//    //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//    Eij = x(N*i+j);        
//    Eijm1 = x(N*i+j-1);
//    Eip1j = x(N*(i-1)+j);
//    Eim1j = x(N*(i-1)+j);
//    Eip1jm1 = x(N*(i-1)+j-1);
//    Eim1jm1 = x(N*(i-1)+j-1);
//    Tij = x(Nsq+N*i+j);
//    Tijp1 = x(Nsq+N*i+j-1);
//    Tijm1 = x(Nsq+N*i+j-1);
//    Tip1j = x(Nsq+N*(i-1)+j);
//    Tim1j = x(Nsq+N*(i-1)+j);
//    Eijp1 = Eijm1-3.*h*sigma(j*h,i*h,Tij)*Eij;
//    Eip1jp1 = Eip1jm1-3.*h*sigma(j*h,(i+1)*h,Tip1j)*Eip1j;
//    Eim1jp1 = Eim1jm1-3.*h*sigma(j*h,(i-1)*h,Tim1j)*Eim1j;
//    Eijp05 = (Eij+Eijp1)/2.;
//    Eijm05 = (Eij+Eijm1)/2.;
//    Eip05j = (Eij+Eip1j)/2.;
//    Eim05j = (Eij+Eim1j)/2.;
//    Tijp05 = (Tij+Tijp1)/2.;
//    Tijm05 = (Tij+Tijm1)/2.;
//    Tip05j = (Tij+Tip1j)/2.;
//    Tim05j = (Tij+Tim1j)/2.;
//    gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//    gradEijp05 = sqrt(gradEijp05);
//    gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//    gradEijm05 = sqrt(gradEijm05);
//    gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//    gradEip05j = sqrt(gradEip05j);
//    gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//    gradEim05j = sqrt(gradEim05j);
//
//    // computing diffusion
//    fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                        -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                        +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                        -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//    fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                            -D2(Tijm05)*(Tij  -Tijm1)
//                            +D2(Tip05j)*(Tip1j-Tij)
//                            -D2(Tim05j)*(Tij  -Tim1j));
//
//    //adding radiation
//    radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//    fx(i*N+j) += radiation;
//    fx(Nsq+i*N+j) -= radiation;
//    
//    //lower right corner
//    //here we use homogeneous Neumann for T and E/4-dE/dx1/6/sigma=0 for E
//    //yielding that value for Eijp1
//    i=0;
//    j=N-1;
//    //computing |grad_E|, E and T in (i,j+1/2),(i,j-1/2),(i+1/2,j),(i-1/2,j)
//    Eij = x(N*i+j);        
//    Eijm1 = x(N*i+j-1);
//    Eip1j = x(N*(i+1)+j);
//    Eim1j = x(N*(i+1)+j);
//    Eip1jm1 = x(N*(i+1)+j-1);
//    Eim1jm1 = x(N*(i+1)+j-1);
//    Tij = x(Nsq+N*i+j);
//    Tijp1 = x(Nsq+N*i+j-1);
//    Tijm1 = x(Nsq+N*i+j-1);
//    Tip1j = x(Nsq+N*(i+1)+j);
//    Tim1j = x(Nsq+N*(i+1)+j);
//    Eijp1 = Eijm1-3.*h*sigma(j*h,i*h,Tij)*Eij;
//    Eip1jp1 = Eip1jm1-3.*h*sigma(j*h,(i+1)*h,Tip1j)*Eip1j;
//    Eim1jp1 = Eim1jm1-3.*h*sigma(j*h,(i-1)*h,Tim1j)*Eim1j;
//    Eijp05 = (Eij+Eijp1)/2.;
//    Eijm05 = (Eij+Eijm1)/2.;
//    Eip05j = (Eij+Eip1j)/2.;
//    Eim05j = (Eij+Eim1j)/2.;
//    Tijp05 = (Tij+Tijp1)/2.;
//    Tijm05 = (Tij+Tijm1)/2.;
//    Tip05j = (Tij+Tip1j)/2.;
//    Tim05j = (Tij+Tim1j)/2.;
//    gradEijp05 = pow((Eijp1-Eij)/h,2)+ pow((Eip1j+Eip1jp1-Eim1j-Eim1jp1)/4/h,2);
//    gradEijp05 = sqrt(gradEijp05);
//    gradEijm05 = pow((Eij-Eijm1)/h,2)+ pow((Eip1j+Eip1jm1-Eim1j-Eim1jm1)/4/h,2);
//    gradEijm05 = sqrt(gradEijm05);
//    gradEip05j = pow((Eip1jp1+Eijp1-Eip1jm1-Eijm1)/4./h,2)+ pow((Eip1j-Eij)/h,2);
//    gradEip05j = sqrt(gradEip05j);
//    gradEim05j = pow((Eim1jp1+Eijp1-Eim1jm1-Eijm1)/4./h,2)+ pow((Eij-Eim1j)/h,2);
//    gradEim05j = sqrt(gradEim05j);
//
//    // computing diffusion
//    fx(i*N+j) = 1./h/h*( D1((j+0.5)*h,i*h,gradEijp05,Eijp05,Tijp05)*(Eijp1-Eij) 
//                        -D1((j-0.5)*h,i*h,gradEijm05,Eijm05,Tijm05)*(Eij-Eijm1)
//                        +D1(j*h,(i+0.5)*h,gradEip05j,Eip05j,Tip05j)*(Eip1j-Eij)
//                        -D1(j*h,(i-0.5)*h,gradEim05j,Eim05j,Tim05j)*(Eij-Eim1j));
//    fx(Nsq+i*N+j) = 1./h/h*( D2(Tijp05)*(Tijp1-Tij) 
//                            -D2(Tijm05)*(Tij  -Tijm1)
//                            +D2(Tip05j)*(Tip1j-Tij)
//                            -D2(Tim05j)*(Tij  -Tim1j));
//
//    //adding radiation
//    radiation = sigma(j*h,i*h,Tij)*(pow(Tij,4)-Eij);
//    fx(i*N+j) += radiation;
//    fx(Nsq+i*N+j) -= radiation;
//}
//
//void IntegroDifferentialEquation::f(Real t, Vector& x, Vector& fx)
//{
//    int i;
//    Real xi;
//    
//    fx(0) = N*N*(x(1)-2.*x(0)+1.-sqrt(t)/2.); // diff with dirichlet 1-sqrt(t)/2
//    for(i=1;i<N-1;i++)
//    {
//        fx(i) = N*N*(x(i-1)-2.*x(i)+x(i+1));
//        
//        xi = (i+1.)*h;
//        fx(i) -= sigma*h/2.*( pow(1.-sqrt(t)/2.,4)/pow(1.+xi,2) + pow(x(N-1),4)/pow(2.-xi,2) );
//        for(int j=1;j<N;j++)
//            fx(i) -= sigma*h*pow(x(j-1),4)/pow(1.+abs(xi-j*h),2);
//    }
//    fx(N-1) = N*N*2.*(x(N-2)-x(N-1)); // diff with neumann
//
//
//    i=0;
//    xi = (i+1.)*h;
//    fx(i) -= sigma*h/2.*( pow(1.-sqrt(t)/2.,4)/pow(1.+xi,2) + pow(x(N-1),4)/pow(2.-xi,2) );
//    for(int j=1;j<N;j++)
//        fx(i) -= sigma*h*pow(x(j-1),4)/pow(1.+abs(xi-j*h),2);
//    
//    i=N-1;
//    xi = (i+1.)*h;
//    fx(i) -= sigma*h/2.*( pow(1.-sqrt(t)/2.,4)/pow(1.+xi,2) + pow(x(N-1),4)/pow(2.-xi,2) );
//    for(int j=1;j<N;j++)
//        fx(i) -= sigma*h*pow(x(j-1),4)/pow(1.+abs(xi-j*h),2);
//    
//}