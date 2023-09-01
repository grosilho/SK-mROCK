#include "exact_sol.h"

HyperbolicSolution::HyperbolicSolution()
{
    alpha=1.;
    beta=1.;
    gamma=1.;
    delta=1.;
    lambda=1.;
}


Real HyperbolicSolution::eval(Real x, Real y) const
{   
    Real X = gamma*cosh(lambda*x)+delta*sinh(lambda*x);
    Real Y = alpha*cos(lambda*y)+beta*sin(lambda*y);
           
    return X*Y;
}

Vector HyperbolicSolution::eval_grad(Real x, Real y) const
{   
    Real X = gamma*cosh(lambda*x)+delta*sinh(lambda*x);
    Real Y = alpha*cos(lambda*y)+beta*sin(lambda*y);
    
    Vector grad(2);
    grad(0) = lambda*(gamma*sinh(lambda*x)+delta*cosh(lambda*x))*Y;
    grad(1) = lambda*X*(-alpha*sin(lambda*y)+beta*cos(lambda*y));
    
    return grad;
}

Ue::Ue(Real si_, Real se_)
{
    si=si_;
    se=se_;
}

Real Ue::eval(Real x, Real y) const
{
    Real rsq = x*x+y*y;
    return (si/se)*(8.+rsq/2.)/(3.*rsq)*y;
}

Vector Ue::eval_grad(Real x, Real y) const
{
    Real rsq = x*x+y*y;
    Vector grad(2);
    grad(0) = -(si/se)*16.*x*y/(3.*rsq*rsq);
    grad(1) = -(si/se)*16.*y*y/(3.*rsq*rsq)
              +(si/se)*(16.+rsq)/(6.*rsq);
    return grad;
}

Real Ui::eval(Real x, Real y) const
{
    return -y/2.;
}

Vector Ui::eval_grad(Real x, Real y) const
{
    Vector grad(2);
    grad(0) = 0.;
    grad(1) = -0.5;
    return grad;
}

Some2dFunction::Some2dFunction()
{
    a=1.;
    b=1.;
}

Real Some2dFunction::eval(Real x, Real y) const
{
    return cos(a*M_PI*x)*sin(b*M_PI*y);
}

Vector Some2dFunction::eval_grad(Real x, Real y) const
{
    Vector grad(2);
    grad(0) = -a*M_PI*sin(a*M_PI*x)*sin(b*M_PI*y);
    grad(1) = b*M_PI*cos(a*M_PI*x)*cos(b*M_PI*y);
    return grad;
}