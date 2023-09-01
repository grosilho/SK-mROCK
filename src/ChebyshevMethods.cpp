#include "ChebyshevMethods.h"


void ChebyshevMethods::CoefficientsRKC1(vector<Real>& mu, vector<Real>& nu, 
                       vector<Real>& kappa, vector<Real>& c, vector<Real>& d,
                       unsigned int s, Real eps)
{
    mu.resize(s);
    nu.resize(s);
    kappa.resize(s);
    c.resize(s);
    d.resize(s);
    
    Real w0 = 1.+eps/s/s;
    Real w1 = T(w0,s)/dT(w0,s);
    
    Real bjm2,bjm1,bj;
    
    bjm1 = 1./w0; // = 1/T(w0,1)
    mu[0] = w1/w0;
    c[0] = mu[0];
    d[0] = 0;
    
    Real Tjw0,Tjm1w0,Tjm2w0;
    Tjm2w0=1;
    Tjm1w0=w0;
    
    for(unsigned int j=2;j<=s;j++)
    {
        Tjw0 = 2.*w0*Tjm1w0-Tjm2w0;
        bj = 1./Tjw0;
        mu[j-1] = 2*w1*bj/bjm1;
        nu[j-1] = 2*w0*bj/bjm1;
        if(j>2)
        {
            kappa[j-1] = -bj/bjm2;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2]+kappa[j-1]*c[j-3];
            d[j-1] = 2*mu[j-1]*c[j-2]+nu[j-1]*d[j-2]+kappa[j-1]*d[j-3];
        }
        else
        {
            kappa[j-1]=-bj;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2];
            d[j-1] = 2*mu[j-1]*c[j-2];
        }
        
        bjm2 = bjm1;
        bjm1 = bj;
        
        Tjm2w0 = Tjm1w0;
        Tjm1w0 = Tjw0;
    }
}

Real ChebyshevMethods::ls_RKC1(unsigned int s, Real eps)
{
    Real w0 = 1.+eps/(s*s);
    Real w1 = T(w0,s)/dT(w0,s);
    return 2*w0/w1;
}

void ChebyshevMethods::CoefficientsRKC2(vector<Real>& mu, vector<Real>& nu, 
                       vector<Real>& kappa, vector<Real>& gamma, vector<Real>& c,
                       unsigned int s, Real eps)
{
    mu.resize(s,0.);
    nu.resize(s,0.);
    kappa.resize(s,0.);
    gamma.resize(s,0.);
    c.resize(s,0.);
    
    Real w0 = 1.+eps/s/s;
    Real w1 = dT(w0,s)/ddT(w0,s);
    
    Real bjm2,bjm1,bj;
    Real Tjw0,Tjm1w0,Tjm2w0;
    Real dTjw0,dTjm1w0,dTjm2w0;
    Real ddTjw0,ddTjm1w0,ddTjm2w0;
    Tjm2w0=1.;
    Tjm1w0=w0;
    dTjm2w0=0.;
    dTjm1w0=1.;
    ddTjm2w0=0.;
    ddTjm1w0=0.;
    
    for(unsigned int j=2;j<=s;j++)
    {
        Tjw0 = 2.*w0*Tjm1w0-Tjm2w0;
        dTjw0 = 2.*w0*dTjm1w0-dTjm2w0+2.*Tjm1w0;
        ddTjw0 = 2.*w0*ddTjm1w0-ddTjm2w0+4.*dTjm1w0;
        
        bj = ddTjw0/dTjw0/dTjw0;
        
        if(j==2)
        {
            bjm1=bj;
            bjm2=bj;
            mu[0]=bjm1*w1;
            c[0]=mu[0];
        }
        
        mu[j-1] = 2*w1*bj/bjm1;
        nu[j-1] = 2*w0*bj/bjm1;
        kappa[j-1] = -bj/bjm2;
        gamma[j-1] = -(1.-bjm1*Tjm1w0)*mu[j-1];
        
        if(j>2)
            c[j-1] = mu[j-1]+gamma[j-1]+nu[j-1]*c[j-2]+kappa[j-1]*c[j-3];
        else// j==2
            c[j-1] = mu[j-1]+gamma[j-1]+nu[j-1]*c[j-2];
        
        bjm2 = bjm1;
        bjm1 = bj;
        Tjm2w0 = Tjm1w0;
        Tjm1w0 = Tjw0;
        dTjm2w0 = dTjm1w0;
        dTjm1w0 = dTjw0;
        ddTjm2w0 = ddTjm1w0;
        ddTjm1w0 = ddTjw0;
    }
}

Real ChebyshevMethods::ls_RKC2(unsigned int s, Real eps)
{
    Real w0 = 1.+eps/(s*s);
    Real w1 = dT(w0,s)/ddT(w0,s);
    return 2*w0/w1;
}

void ChebyshevMethods::CoefficientsSKROCK(vector<Real>& mu, vector<Real>& nu, 
                       vector<Real>& kappa, vector<Real>& c,unsigned int s, Real eps)
{
    mu.resize(s);
    nu.resize(s);
    kappa.resize(s);
    c.resize(s);
    
    Real w0 = 1.+eps/s/s;
    Real w1 = T(w0,s)/dT(w0,s);
    
    Real bjm2,bjm1,bj;
    
    bjm1 = 1./w0; // = 1/T(w0,1)
    mu[0] = w1/w0;
    nu[0] = s*w1/2./w0;
    kappa[0] = s*w1/w0;
    c[0] = mu[0];
    
    Real Tjw0,Tjm1w0,Tjm2w0;
    Tjm2w0=1;
    Tjm1w0=w0;
    
    for(unsigned int j=2;j<=s;j++)
    {
        Tjw0 = 2.*w0*Tjm1w0-Tjm2w0;
        bj = 1./Tjw0;
        mu[j-1] = 2*w1*bj/bjm1;
        nu[j-1] = 2*w0*bj/bjm1;
        if(j>2)
        {
            kappa[j-1] = -bj/bjm2;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2]+kappa[j-1]*c[j-3];
        }
        else
        {
            kappa[j-1]=-bj;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2];
        }
        
        bjm2 = bjm1;
        bjm1 = bj;
        
        Tjm2w0 = Tjm1w0;
        Tjm1w0 = Tjw0;
    }
}

Real ChebyshevMethods::ls_SKROCK(unsigned int s, Real eps)
{
    return ls_RKC1(s,eps);
}

void ChebyshevMethods::CoefficientsRKU1(vector<Real>& mu, vector<Real>& nu, 
                       vector<Real>& kappa, vector<Real>& c, vector<Real>& d,
                       unsigned int s, Real eps)
{
    mu.resize(s);
    nu.resize(s);
    kappa.resize(s);
    c.resize(s);
    d.resize(s);
    
    Real w0 = 1.+3.*eps/(s*(s+1.)*(s+2.));
    Real w1 = U(w0,s)/dU(w0,s);
    
    Real bjm2,bjm1,bj;
    
    bjm1 = 1./(2.*w0); // = 1/U(w0,1)
    mu[0] = w1/w0;
    c[0] = mu[0];
    d[0] = 0;
    
    Real Ujw0,Ujm1w0,Ujm2w0;
    Ujm2w0=1.;
    Ujm1w0=2.*w0;
    
    for(unsigned int j=2;j<=s;j++)
    {
        Ujw0 = 2.*w0*Ujm1w0-Ujm2w0;
        bj = 1./Ujw0;
        mu[j-1] = 2*w1*bj/bjm1;
        nu[j-1] = 2*w0*bj/bjm1;
        if(j>2)
        {
            kappa[j-1] = -bj/bjm2;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2]+kappa[j-1]*c[j-3];
            d[j-1] = 2*mu[j-1]*c[j-2]+nu[j-1]*d[j-2]+kappa[j-1]*d[j-3];
        }
        else
        {
            kappa[j-1]=-bj;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2];
            d[j-1] = 2*mu[j-1]*c[j-2];
        }
        
        bjm2 = bjm1;
        bjm1 = bj;
        
        Ujm2w0 = Ujm1w0;
        Ujm1w0 = Ujw0;
    }
}

Real ChebyshevMethods::ls_RKU1(unsigned int s, Real eps)
{
    Real w0 = 1.+3.*eps/(s*(s+1.)*(s+2.));
    Real w1 = U(w0,s)/dU(w0,s);
    return 2.*w0/w1;
}

void ChebyshevMethods::CoefficientsRKU2(vector<Real>& mu, vector<Real>& nu, 
                       vector<Real>& kappa, vector<Real>& gamma, vector<Real>& c,
                       unsigned int s, Real eps)
{
    mu.resize(s,0.);
    nu.resize(s,0.);
    kappa.resize(s,0.);
    gamma.resize(s,0.);
    c.resize(s,0.);
    
    Real w0 = 1.+3.*eps/(s*(s+1.)*(s+2.));
    Real w1 = dU(w0,s)/ddU(w0,s);
    
    Real bjm2,bjm1,bj;
    Real Ujw0,Ujm1w0,Ujm2w0;
    Real dUjw0,dUjm1w0,dUjm2w0;
    Real ddUjw0,ddUjm1w0,ddUjm2w0;
    Ujm2w0=1.;
    Ujm1w0=2*w0;
    dUjm2w0=0.;
    dUjm1w0=2.;
    ddUjm2w0=0.;
    ddUjm1w0=0.;
    
    for(unsigned int j=2;j<=s;j++)
    {
        Ujw0 = 2.*w0*Ujm1w0-Ujm2w0;
        dUjw0 = 2.*w0*dUjm1w0-dUjm2w0+2.*Ujm1w0;
        ddUjw0 = 2.*w0*ddUjm1w0-ddUjm2w0+4.*dUjm1w0;
        
        bj = ddUjw0/dUjw0/dUjw0;
        
        if(j==2)
        {
            bjm1=bj;
            bjm2=bj;
            mu[0]=2.*bjm1*w1;
            c[0]=mu[0];
        }
        
        mu[j-1] = 2*w1*bj/bjm1;
        nu[j-1] = 2*w0*bj/bjm1;
        kappa[j-1] = -bj/bjm2;
        gamma[j-1] = -(1.-bjm1*Ujm1w0)*mu[j-1];
        
        if(j>2)
            c[j-1] = mu[j-1]+gamma[j-1]+nu[j-1]*c[j-2]+kappa[j-1]*c[j-3];
        else// j==2
            c[j-1] = mu[j-1]+gamma[j-1]+nu[j-1]*c[j-2];
        
        bjm2 = bjm1;
        bjm1 = bj;
        Ujm2w0 = Ujm1w0;
        Ujm1w0 = Ujw0;
        dUjm2w0 = dUjm1w0;
        dUjm1w0 = dUjw0;
        ddUjm2w0 = ddUjm1w0;
        ddUjm1w0 = ddUjw0;
    }
}

Real ChebyshevMethods::ls_RKU2(unsigned int s, Real eps)
{
    Real w0 = 1.+3.*eps/(s*(s+1.)*(s+2.));
    Real w1 = dU(w0,s)/ddU(w0,s);
    return 2.*w0/w1;
}

Real ChebyshevMethods::T(Real x, unsigned int s)
{
    if(s==0)
        return 1.;
    else if(s==1)
        return x;
    else
    {
        Real Tjm2 = 1.;
        Real Tjm1 = x;
        Real Tj;
        for(unsigned int j=2;j<=s;j++)
        {
            Tj = 2*x*Tjm1-Tjm2;
            Tjm2=Tjm1;
            Tjm1=Tj;
        }
        return Tj;
    }
}

Real ChebyshevMethods::U(Real x, unsigned int s)
{
    if(s==0)
        return 1.;
    else if(s==1)
        return 2.*x;
    else
    {
        Real Ujm2 = 1.;
        Real Ujm1 = 2.*x;
        Real Uj;
        for(unsigned int j=2;j<=s;j++)
        {
            Uj = 2*x*Ujm1-Ujm2;
            Ujm2=Ujm1;
            Ujm1=Uj;
        }
        return Uj;
    }
}

Real ChebyshevMethods::dT(Real x, unsigned int s)
{
    if(s==0)
        return 0.;
    else if(s==1)
        return 1.;
    else
        return s*U(x,s-1);
}

Real ChebyshevMethods::dU(Real x, unsigned int s)
{
    if(s==0)
        return 0.;
    else if(s==1)
        return 2;
    else
        return ddT(x,s+1)/(s+1);
}

Real ChebyshevMethods::ddT(Real x, unsigned int s)
{
    if(s==0)
        return 0;
    else if(s==1)
        return 0;
    else
    {
        if(x==1)
            return s*s*(s*s-1)/3.;
        else if(x==-1)
            return pow(-1,s)*s*s*(s*s-1)/3.;
        else
            return s/(x*x-1)*((s+1)*T(x,s)-U(x,s));
    }
}

Real ChebyshevMethods::ddU(Real x, unsigned int s)
{
    if(s==0)
        return 0.;
    else if(s==1)
        return 0.;
    else
    {
        Real Ujm2 = 1.;
        Real Ujm1 = 2.*x;
        Real dUjm2 = 0.;
        Real dUjm1 = 2.;
        Real ddUjm2 = 0.;
        Real ddUjm1 = 0.;
        Real Uj,dUj,ddUj;
        for(unsigned int j=2;j<=s;j++)
        {
            Uj = 2*x*Ujm1-Ujm2;
            dUj = 2*(Ujm1+x*dUjm1)-dUjm2;
            ddUj = 2*(2.*dUjm1+x*ddUjm1)-ddUjm2;
            Ujm2=Ujm1;
            Ujm1=Uj;
            dUjm2=dUjm1;
            dUjm1=dUj;
            ddUjm2=ddUjm1;
            ddUjm1=ddUj;
        }
        return ddUj;
    }
}