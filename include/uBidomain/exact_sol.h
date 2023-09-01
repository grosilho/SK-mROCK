#ifndef EXACT_SOL_H
#define EXACT_SOL_H

#include <iostream>
#include <cmath>
#include<Eigen/Dense>
#include <vector>
#include <iomanip> 


typedef double Real;
typedef Eigen::VectorXd Vector;
typedef Eigen::RowVector2d Point;
typedef Eigen::MatrixXd Matrix;
using namespace std;

class HyperbolicSolution
{
public:
    HyperbolicSolution();
    
    Real eval(Real x, Real y) const;
    Vector eval_grad(Real x, Real y) const;
    
protected:
    Real alpha;
    Real beta;
    Real gamma;
    Real delta;
    Real lambda;
};

class Ue
{
public:
    Ue(Real si_, Real se_);
    Real eval(Real x, Real y) const;
    Vector eval_grad(Real x, Real y) const;
    
protected:
    Real si,se;
};

class Ui
{
public:
    Real eval(Real x, Real y) const;
    Vector eval_grad(Real x, Real y) const;
};

class Some2dFunction
{
public:
    Some2dFunction();
    Real eval(Real x, Real y) const;
    Vector eval_grad(Real x, Real y) const;
    
private:
    Real a;
    Real b;
};


#endif /* EXACT_SOL_H */

