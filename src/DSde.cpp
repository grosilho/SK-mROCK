#include <fstream>

#include "DSde.h"

DSde::DSde()
:Ode()
{
    noise=GENERAL;
    Wsize=0;
}

DSde::~DSde()
{
}

int DSde::brownian_size()
{
    return Wsize;
}

NoiseType DSde::get_noise_type() const
{
    return noise;
}

void DSde::g(Real t, Vector& x, Vector& G, int r)
{
    cout<<"ERROR: using undefined function g(Real, Vector&, Vector&, int)"<<endl;
}

void DSde::g(Real t, Vector& x, Vector& G)
{
    cout<<"ERROR: using undefined function g(Real, Vector&, Vector&)"<<endl;
}

void DSde::g(Real t, Vector& x, Matrix& G)
{
    cout<<"ERROR: using undefined function g(Real, Vector&, Matrix&)"<<endl;
}

MultirateDSde::MultirateDSde()
: Ode(), MultirateOde(), DSde()
{

}

MultirateDSde::~MultirateDSde()
{

}