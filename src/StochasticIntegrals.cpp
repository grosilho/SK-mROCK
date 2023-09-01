#include <chrono>
#include "StochasticIntegrals.h"
#include <algorithm>  // generate
#include <iterator>   // begin, end, and ostream_iterator
#include <functional> // bind
#ifdef _OPENMP
#include <omp.h>
#endif

HistoryStochasticIntegrals::HistoryStochasticIntegrals(int size, 
                            bool continuous_, bool doubleint_,
                            NoiseType noisetype_, int seed)
:m(size), continuous(continuous_), doubleintegral(doubleint_), noisetype(noisetype_)
{   
    if(seed==-1)
    {
#ifdef _OPENMP
    seed = chrono::system_clock::now().time_since_epoch().count()
                    +omp_get_thread_num();
#else
    seed = chrono::system_clock::now().time_since_epoch().count();
#endif
    }
    
    gen.seed(seed);
}

HistoryStochasticIntegrals::~HistoryStochasticIntegrals()
{
    
}

void HistoryStochasticIntegrals::sample_integrals(Real tend, unsigned int npoints)
{
    Real h = tend/npoints;
    I.resize(npoints);
    for(unsigned int i=0;i<npoints;i++)
        I[i] = StochasticIntegral(m,continuous,doubleintegral,noisetype,h,gen);
    
}

void HistoryStochasticIntegrals::sample_integrals(Real tend, Real h)
{
    unsigned int npoints = ceil(tend/h);
    I.resize(npoints);
    Real t=0.;
    
    for(unsigned int i=0;i<npoints;i++)
    {
        if(t+h>tend)
            h=tend-t;
        
        t+=h;
        I[i] = StochasticIntegral(m,continuous,doubleintegral,noisetype,h,gen);
    }
    
}

void HistoryStochasticIntegrals::coarse()
{
    unsigned int npoints_new = I.size()/2;
    if(npoints_new*2!=I.size())
    {
        cout<<"ERROR: cannot coarse BrownianMotion anymore"<<endl;
        return;
    }
    
    for(unsigned int i=0;i<npoints_new;i++)
        I[i] = I[2*i]+I[2*i+1];
    I.resize(npoints_new);
}

void HistoryStochasticIntegrals::write_paths(string filename)
{
    Vector W(m);
    W *= 0.;
    Real t=0.;
    ofstream file(filename+string(".csv"),ios::out);
    
    for(unsigned int i=0;i<I.size();i++)
    {
        W += I[i].get_Ir();
        t += I[i].get_h();
        for(unsigned int j=0;j<m;j++)
            file<<W(j)<<", ";
        file<<t<<endl;     
    }
    file.close();
}

int HistoryStochasticIntegrals::integral_size()
{
    return m;
}

int HistoryStochasticIntegrals::length()
{
    return I.size();
}

StochasticIntegral& HistoryStochasticIntegrals::operator [](const unsigned int i)
{
    return I[i];
}

// -----------------------------------------------------
StochasticIntegral::StochasticIntegral()
{
}

StochasticIntegral::StochasticIntegral(unsigned int m_, bool continuous_, 
                    bool doubleintegral_, NoiseType noisetype_, Real h_)
:m(m_),Ir(m_),h(h_),
 continuous(continuous_),doubleintegral(doubleintegral_),noisetype(noisetype_)
{
    if(doubleintegral)
        Ipq = Matrix::Zero(m+1,m+1);
}

StochasticIntegral::StochasticIntegral(unsigned int m_, bool continuous_,
                    bool doubleintegral_, NoiseType noisetype_, Real h_, default_random_engine& gen)
:m(m_),Ir(m_),h(h_),
 continuous(continuous_),doubleintegral(doubleintegral_),noisetype(noisetype_),
 normal(0.0,1.0),xid(1,6)
{
    if(doubleintegral)
        Ipq = Matrix::Zero(m+1,m+1); // integrals \int_0^h\int_0^t dWp(s) dWq(t). For p=0 then dWp(s)=ds.
    
    if(continuous)
        sample_continuous(gen);
    else
        sample_discrete(gen);
}

void StochasticIntegral::sample_continuous(default_random_engine& gen)
{
    Real sqrh = sqrt(h);
    
    if(!doubleintegral)
    {
        for(int i=0;i<m;i++)
            Ir(i) = sqrh*normal(gen);
    }
    else if(noisetype==DIAGONAL || noisetype==COMMUTATIVE)
    {
        Real xi;
        for(int i=0;i<m;i++)
        {
            xi = normal(gen);
            
            Ir(i) = sqrh*xi;
            Ipq(i+1,i+1) = h*(xi*xi-1.)/2.;
        }
        
        if(noisetype==COMMUTATIVE)// we suoppose that Ipq and Iqp will be used summed together, hence we give to each one the same contribution
            for(int i=1;i<=m;i++)
            for(int j=i+1;j<=m;j++)
            {
                Ipq(i,j) = Ir[i-1]*Ir[j-1]/2.;
                Ipq(j,i) = Ir[i-1]*Ir[j-1]/2.;
            }
    }
    else if(noisetype==GENERAL)
    {
        const unsigned int N = floor(2./h);
        Real xi[m];
        Real mu[m];
        Real nu[m][N];
        Real zeta[m][N];
        for(unsigned int j=0;j<m;j++)
        {
            xi[j] = normal(gen);
            mu[j] = normal(gen);
            for(unsigned int r=0;r<N;r++)
            {
                nu[j][r] = normal(gen);
                zeta[j][r] = normal(gen);
            }
        }
        
        const Real Pi = 3.14159265358979323846;
        const Real sqr2 = sqrt(2.);
        Real sqrRhoN=0.;
        for(int r=1;r<=N;r++)
            sqrRhoN += (1./r)/r;
        sqrRhoN = 1./12.-1./2./Pi/Pi*sqrRhoN;
        sqrRhoN = sqrt(sqrRhoN);

        for(unsigned int i=0;i<m;i++)
        {
            Ir(i) = sqrh*xi[i];
            Ipq(i+1,i+1) = h*(xi[i]*xi[i]-1.)/2.;
            
            for(int j=i+1;j<m;j++)
            {
                Ipq(i+1,j+1) = 0.;
                for(unsigned int r=0;r<N;r++)
                    Ipq(i+1,j+1) += (zeta[i][r]*(sqr2*xi[j]+nu[j][r])-zeta[j][r]*(sqr2*xi[i]+nu[i][r]))/(r+1.);
                Ipq(i+1,j+1) *= h/2./Pi;
                Ipq(i+1,j+1) += h*(xi[i]*xi[j]/2.+sqrRhoN*(mu[i]*xi[j]-mu[j]*xi[i]));
                Ipq(j+1,i+1) = h*xi[i]*xi[j]-Ipq(i+1,j+1);
            }
        }
        
    }
}

void StochasticIntegral::sample_discrete(default_random_engine& gen)
{
    Real sqrh = sqrt(h);
    Real sqr3 = sqrt(3.);
    int tmp;
    
    if(!doubleintegral)
    {
        Real xi;
        for(int i=0;i<m;i++)
        {
            tmp = xid(gen);
            xi = (tmp<=4.) ? 0. : ((tmp==5.) ? sqr3 : -sqr3) ;

            Ir(i) = sqrh*xi;
        }
    }
    else if(noisetype==DIAGONAL || noisetype==COMMUTATIVE) // I think that with this we are ok up to weak order 2
    {
        Real xi;
        Ipq(0,0) = h*h/2.;
        for(int i=0;i<m;i++)
        {
            tmp = xid(gen);
            xi = (tmp<=4.) ? 0. : ((tmp==5.) ? sqr3 : -sqr3) ;
            
            Ir(i) = sqrh*xi;
            Ipq(i+1,i+1) = h*(xi*xi-1.)/2.;
            Ipq(0,i+1) = h*sqrh*xi/2.;
            Ipq(i+1,0) = h*sqrh*xi/2.;
        }
        
        if(noisetype==COMMUTATIVE)
            for(int i=1;i<=m;i++)
            for(int j=i+1;j<=m;j++)
            {
                Ipq(i,j) = Ir[i-1]*Ir[j-1]/2.;
                Ipq(j,i) = Ir[i-1]*Ir[j-1]/2.;
            }
        
    }
    else if(noisetype==GENERAL)
    {
        Real xi[m];
        Real chi[m];
        
        for(unsigned int j=0;j<m;j++)
        {
            tmp = xid(gen);
            xi[j] = (tmp<=4.) ? 0. : ((tmp==5.) ? sqr3 : -sqr3) ;
            chi[j] = (xid(gen)<=3) ? 1:-1;
        }
        
        for(int i=0;i<m;i++)
        {
            Ir(i) = sqrh*xi[i];
            Ipq(i+1,i+1) = h*(xi[i]*xi[i]-1.)/2.;
            Ipq(0,i+1) = h*sqrh*xi[i]/2.;
            Ipq(i+1,0) = h*sqrh*xi[i]/2.;
            
            for(int j=i+1;j<m;j++)
            {
                Ipq(i+1,j+1) = h*(xi[i]*xi[j]+chi[i]*chi[j])/2.; 
                Ipq(j+1,i+1) = h*(xi[i]*xi[j]-chi[i]*chi[j])/2.;
            }
        }
    }
}

Vector& StochasticIntegral::get_Ir()
{
    return Ir;
}

Matrix& StochasticIntegral::get_Ipq()
{
    return Ipq;
}

Real StochasticIntegral::get_h()
{
    return h;
}

StochasticIntegral StochasticIntegral::operator +(const StochasticIntegral& si)
{
    if(si.m!=this->m)
        cout<<" ERROR: adding two StochasticIntegrals of different sizes."<<endl;
    
    StochasticIntegral newsi(this->m,this->continuous,
                             this->doubleintegral,this->noisetype,this->h+si.h);
    
    newsi.Ir = this->Ir + si.Ir;
    
    if(doubleintegral)
    {
        newsi.Ipq(0,0) = 0.5*newsi.h*newsi.h;
        for(int j=1;j<=m;j++)
        {
            newsi.Ipq(0,j) = this->Ipq(0,j)+si.Ipq(0,j)+this->h*si.Ir(j-1);
            newsi.Ipq(j,0) = this->Ipq(j,0)+si.Ipq(j,0)+this->Ir(j-1)*si.h;
            
            for(int i=1;i<=m;i++)
                newsi.Ipq(i,j) = this->Ipq(i,j)+si.Ipq(i,j)+this->Ir(i-1)*si.Ir(j-1);
        }
    }
    
    return newsi;
}