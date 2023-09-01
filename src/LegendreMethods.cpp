#include "LegendreMethods.h"

void LegendreMethods::CoefficientsRKL1(vector<Real>& mu, vector<Real>& nu, 
                       vector<Real>& kappa, vector<Real>& c, vector<Real>& d,
                       unsigned int s, Real eps)
{
    mu.resize(s);
    nu.resize(s);
    kappa.resize(s);
    c.resize(s);
    d.resize(s);
    
    Real w0 = 1.+2.*eps/(s*s+s);
    Real w1 = L(w0,s)/dL(w0,s);
    
    Real bjm2,bjm1,bj;
    
    bjm1 = 1./w0; // = 1/T(w0,1)
    mu[0] = w1/w0;
    c[0] = mu[0];
    d[0] = 0;
    
    Real Ljw0,Ljm1w0,Ljm2w0;
    Ljm2w0=1.;
    Ljm1w0=w0;
    
    for(unsigned int j=2;j<=s;j++)
    {
        Ljw0 = ((2.*j-1.)*w0*Ljm1w0-(j-1.)*Ljm2w0)/j;
        bj = 1./Ljw0;
        mu[j-1] = ((2.*j-1.)/j)*w1*bj/bjm1;
        nu[j-1] = ((2.*j-1.)/j)*w0*bj/bjm1;
        if(j>2)
        {
            kappa[j-1] = -((j-1.)/j)*bj/bjm2;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2]+kappa[j-1]*c[j-3];
            d[j-1] = 2*mu[j-1]*c[j-2]+nu[j-1]*d[j-2]+kappa[j-1]*d[j-3];
        }
        else
        {
            kappa[j-1]= -((j-1.)/j)*bj;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2];
            d[j-1] = 2*mu[j-1]*c[j-2];
        }
        
        bjm2 = bjm1;
        bjm1 = bj;
        
        Ljm2w0 = Ljm1w0;
        Ljm1w0 = Ljw0;
    }
}

void LegendreMethods::CoefficientsRKL2(vector<Real>& mu, vector<Real>& nu, 
                       vector<Real>& kappa, vector<Real>& gamma, vector<Real>& c,
                       unsigned int s, Real eps)
{   
    mu.resize(s,0.);
    nu.resize(s,0.);
    kappa.resize(s,0.);
    gamma.resize(s,0.);
    c.resize(s,0.);
    
    Real w0 = 1.+2.*eps/s/(s+1.);
    Real w1 = dL(w0,s)/ddL(w0,s);
        
    Real bjm2,bjm1,bj;
    Real Ljw0,Ljm1w0,Ljm2w0;
    Real dLjw0,dLjm1w0,dLjm2w0;
    Real ddLjw0,ddLjm1w0,ddLjm2w0;
    
    Ljm2w0=1.;
    Ljm1w0=w0;
    dLjm2w0=0.;
    dLjm1w0=1.;
    ddLjm2w0=0.;
    ddLjm1w0=0.;
    
    for(unsigned int j=2;j<=s;j++)
    {
        Ljw0 = ((2*j-1.)*w0*Ljm1w0-(j-1.)*Ljm2w0)/j;
        dLjw0 = ((2*j-1.)*(w0*dLjm1w0+Ljm1w0)-(j-1.)*dLjm2w0)/j;
        ddLjw0 = (2.*(2*j-1.)*dLjm1w0+(2*j-1.)*w0*ddLjm1w0-(j-1.)*ddLjm2w0)/j;
        
        bj = ddLjw0/(dLjw0*dLjw0);
        
        if(j==2)
        {
            bjm1=bj;
            bjm2=bj;
            mu[0]=bjm1*w1;
            c[0]=mu[0];
        }
        
        mu[j-1] = ((2.*j-1.)/j)*w1*bj/bjm1;
        nu[j-1] = ((2.*j-1.)/j)*w0*bj/bjm1;
        kappa[j-1] = -((j-1.)/j)*bj/bjm2;
        gamma[j-1] = -(1.-bjm1*Ljm1w0)*mu[j-1];
        
        if(j>2)
            c[j-1] = mu[j-1]+gamma[j-1]+nu[j-1]*c[j-2]+kappa[j-1]*c[j-3];
        else// j==2
            c[j-1] = mu[j-1]+gamma[j-1]+nu[j-1]*c[j-2];
        
        bjm2 = bjm1;
        bjm1 = bj;
        Ljm2w0 = Ljm1w0;
        Ljm1w0 = Ljw0;
        dLjm2w0 = dLjm1w0;
        dLjm1w0 = dLjw0;
        ddLjm2w0 = ddLjm1w0;
        ddLjm1w0 = ddLjw0;
    }
}

Real LegendreMethods::ls_RKL1(unsigned int s, Real eps)
{
    Real w0 = 1.+2.*eps/(s*s+s);
    Real w1 = L(w0,s)/dL(w0,s);
    return 2*w0/w1;
}

Real LegendreMethods::ls_RKL2(unsigned int s, Real eps)
{
    Real w0 = 1.+2.*eps/s/(s+1.);
    Real w1 = dL(w0,s)/ddL(w0,s);
    return 2*w0/w1;
}

Real LegendreMethods::L(Real x, unsigned int s)
{
    if(s==0)
        return 1.;
    else if(s==1)
        return x;
    else
    {
        Real Pjm2 = 1.;
        Real Pjm1 = x;
        Real Pj;
        for(unsigned int j=2;j<=s;j++)
        {
            Pj = ((2*j-1.)*x*Pjm1-(j-1.)*Pjm2)/j;
            Pjm2=Pjm1;
            Pjm1=Pj;
        }
        return Pj;
    }
}

Real LegendreMethods::dL(Real x, unsigned int s)
{
    if(s==0)
        return 0.;
    else if(s==1)
        return 1.;
    else
    {
        Real Pjm2 = 1.;
        Real Pjm1 = x;
        Real dPjm2 = 0.;
        Real dPjm1 = 1.;
        Real Pj, dPj;
        for(unsigned int j=2;j<=s;j++)
        {
            Pj = ((2*j-1.)*x*Pjm1-(j-1.)*Pjm2)/j;
            dPj = ((2*j-1.)*(x*dPjm1+Pjm1)-(j-1.)*dPjm2)/j;
            
            Pjm2=Pjm1;
            Pjm1=Pj;
            dPjm2=dPjm1;
            dPjm1=dPj;
        }
        return dPj;
    }
}

Real LegendreMethods::ddL(Real x, unsigned int s)
{
    if(s==0)
        return 0.;
    else if(s==1)
        return 0.;
    else
    {
        Real Pjm2 = 1.;
        Real Pjm1 = x;
        Real dPjm2 = 0.;
        Real dPjm1 = 1.;
        Real ddPjm2 = 0.;
        Real ddPjm1 = 0.;
        Real Pj, dPj, ddPj;
        
        for(unsigned int j=2;j<=s;j++)
        {
            Pj = ((2*j-1.)*x*Pjm1-(j-1.)*Pjm2)/j;
            dPj = ((2*j-1.)*(x*dPjm1+Pjm1)-(j-1.)*dPjm2)/j;
            ddPj = (2.*(2*j-1.)*dPjm1+(2*j-1.)*x*ddPjm1-(j-1.)*ddPjm2)/j;
            
            Pjm2=Pjm1;
            Pjm1=Pj;
            dPjm2=dPjm1;
            dPjm1=dPj;
            ddPjm2=ddPjm1;
            ddPjm1=ddPj;
        }
        return ddPj;
    }
}