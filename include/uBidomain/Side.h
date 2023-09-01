#ifndef SIDE_H
#define SIDE_H

#include <utility>
#include <cmath>
#include <string>
#include <Eigen/Dense>
#include <iostream>

using namespace std;

typedef double Real;

typedef Eigen::MatrixXd Matrix;
typedef Eigen::RowVector2d Point;
typedef Eigen::VectorXd Vector;
typedef Eigen::ArrayXd Array;

namespace uBidomain
{

class Side
{
public:
    Side();
    Side(const vector<Side>& sides);
    
    void resize_vectors(unsigned int n);
    void concatenate_and_set(const vector<Side>& sides);
    
    const Vector& get_x() const;
    const Matrix& get_G() const;
    const Matrix& get_dG() const;
    const Matrix& get_ddG() const;
    const Matrix& get_N() const;
    const Vector& get_norm_dG() const;
    const Vector& get_norm_ddG() const;
    
    const Vector& get_kappa() const;
    
    Point get_G(unsigned int i) const;    
    Point get_N(unsigned int i) const;
    Real get_norm_dG(unsigned int i) const;
   
    unsigned int get_n() const;
    Real get_dx() const;
    Real get_h() const;
    Real get_length() const;
    
    Point begin() const;
    Point end() const;
    
    ostream& write(ostream& os, string name="");
    
    Side reverse_normal();
    Side operator-() const;
    Side operator+(const Point& p) const;
    Side operator*(const Real r) const;
    
protected:
    Vector x;
    Matrix G;
    Matrix dG;
    Matrix ddG;
    Matrix N;
    Vector norm_dG;
    Vector norm_ddG;
    
    Point p1;
    Point p2;
    
    unsigned int n;
    Real dx;
    Real length;
    Real h;
    
    Vector kappa;
};

class Segment: public Side
{
public:
    Segment(Point p1_, Point p2_, Real h_, Real sign_N, Real kappa_);
};

class Arc: public Side
{
public: 
    Arc(Point c_, Real r_, Real alpha_, Real beta_, Real h_, Real sign_N, Real kappa_);
    
protected:
    Point c;
    Real r;
    Real alpha;
    Real beta;
};

class EllipseArc: public Side
{
public: 
    EllipseArc(Point c_, Real r1_, Real r2_, Real alpha_, Real beta_, Real h_, Real sign_N, Real kappa_);
    
protected:
    Point c;
    Real r1;
    Real r2;
    Real alpha;
    Real beta;
};

class Wave1: public Side
{
public:
    Wave1(Point p1_, Point p2_, unsigned k_, Real a_, Real h_, Real sign_N, Real kappa_);
    
protected:
    unsigned k;
    Real a;
};

class Wave2: public Side
{
public:
    Wave2(Point p1_, Point p2_, unsigned k_, Real a_, Real h_, Real sign_N, Real kappa_);
    
protected:
    unsigned k;
    Real a;
};

class SquaredWave1: public Side
{
public:
    SquaredWave1(Point p1_, Point p2_, unsigned k_, Real a_, Real h_, Real sign_N, Real kappa_dir, Real kappa_perp);
    
protected:
    unsigned k;
    Real a;
};

class SquaredWave2: public Side
{
public:
    SquaredWave2(Point p1_, Point p2_, unsigned k_, Real a_, Real h_, Real sign_N, Real kappa_dir, Real kappa_perp);
    
protected:
    unsigned k;
    Real a;
};

} // end namespace uBidomain

ostream& operator<<(ostream& os, const Matrix& M);
void write_matrix(ostream& os, const Matrix& M, string name);
ostream& operator<<(ostream& os, const uBidomain::Side& s);
uBidomain::Side operator*(const Real r,const uBidomain::Side& s);
uBidomain::Side operator+(const Point p,const uBidomain::Side& s);
uBidomain::Side operator-(const Point p,const uBidomain::Side& s);

#endif /* SIDE_H */
