#ifndef STOCHASTICINTEGRAL_H
#define	STOCHASTICINTEGRAL_H

#include "MainHeader.h"

enum NoiseType {DIAGONAL, COMMUTATIVE, GENERAL};

class StochasticIntegral;

class HistoryStochasticIntegrals
{
public:
    HistoryStochasticIntegrals(int size, bool continuous_, bool doubleint_, 
                               NoiseType noisetype_, int seed=-1);
    virtual ~HistoryStochasticIntegrals();
    
    virtual void sample_integrals(Real tend, unsigned int npoints);
    virtual void sample_integrals(Real tend, Real h);
    virtual void coarse();

    int integral_size();
    int length();
    void write_paths(string filename);
    
    StochasticIntegral& operator[](const unsigned int i);
    
protected:
    int m;     // number of independent brownian motions W^i

    vector<StochasticIntegral> I; // history of integrals
    
    bool continuous;
    const bool doubleintegral;
    NoiseType noisetype;
    
    default_random_engine gen;
};

class StochasticIntegral
{
public:
    StochasticIntegral();
    StochasticIntegral(unsigned int m_, bool continuous_,
                       bool doubleintegral_, NoiseType noisetype_,
                       Real h_);
    StochasticIntegral(unsigned int m_, bool continuous_,
                       bool doubleintegral_, NoiseType noisetype_, 
                       Real h_, default_random_engine& gen);
    
    Vector& get_Ir();
    Matrix& get_Ipq();
    Real get_h();
    
    StochasticIntegral operator+(const StochasticIntegral& si);
          
private:
    void sample_continuous(default_random_engine& gen);
    void sample_discrete(default_random_engine& gen);
    
    unsigned int m;
    Real h;
    Vector Ir;
    Matrix Ipq;
    bool continuous;
    bool doubleintegral;
    NoiseType noisetype;
    normal_distribution<Real> normal;
    uniform_int_distribution<int> xid;
};

#endif	/* STOCHASTICINTEGRAL_H */

