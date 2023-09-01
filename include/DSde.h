#ifndef SDE_H
#define	SDE_H

#include "Ode.h"
#include "StochasticIntegrals.h"


class DSde: public virtual Ode
{
public:    
    DSde();
    virtual ~DSde();
    
    virtual void g(Real t, Vector& x, Vector& G);//only for diagonal noise
    virtual void g(Real t, Vector& x, Matrix& G);//only for general noise
    virtual void g(Real t, Vector& x, Vector& G, int r);//only for general noise
    
    int brownian_size();
    NoiseType get_noise_type() const;
    
    virtual Real phi(const Vector& X) =0;
      
protected:
    int Wsize;
    NoiseType noise;
    
};

#endif	/* SDE_H */

// For diagonal noise
// g(...,Vector) returns a vector where each element 
// is g^i(x^i) for i=1,...,neqn=Wsize

//For general noise
// g(...,Vector,int r) returns the rth diffusion g^r 
// g(...,Matrix) return a matrix with Wsize columns, where each column is g^r(x))
