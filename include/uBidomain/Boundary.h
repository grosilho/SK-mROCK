#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <memory>
#include <vector>
#include <array>
#include <valarray>
#include <map>
#include <set>

#include "Side.h"
#include "FSeries.h"

namespace uBidomain
{
    
class Boundary
{
public:
    Boundary();
    Boundary(vector<Side> sides);
    Boundary(string filename);
    Boundary(ifstream& in);
    
    static Boundary get_circle(Point c, Real r, Real dG, Real sign_N, Real kappa);
    static Boundary get_square(Point c, Real s, Real dG, Real sign_N, Real kappa);
    static Boundary get_rectangle(Point p1, Point p2, Real dG, Real sign_N, Real kappa);
    static Boundary get_ellipse(Point c, Real r1, Real r2, Real dG, Real sign_N, Real kappa);
    
    const unsigned get_n_pts() const;
    const Vector& get_t() const;
    const Vector& get_dt() const;
    Real get_h() const;
    Real get_length() const;
    const Matrix& get_G() const;
    const Point get_G(unsigned i) const;
    const Matrix& get_dG() const;
    const Matrix& get_ddG() const;
    const Matrix& get_N() const;
    const Point get_N(unsigned i) const;
    const Vector& get_norm_dG() const;
    const Vector& get_norm_ddG() const;
    
    const Vector& get_measure() const;
    const Real get_measure(unsigned i) const;
    
    vector<array<unsigned,3>> get_sides_dofs() const;
    vector<array<unsigned,2>> get_corners() const;
    
    Vector interpolate(const Vector& v, const Boundary& bnd) const;
    Real integrate(const Vector& v) const;
    Real L2norm(const Vector& v) const;
    
    ostream& to_matlab(ostream& os, string name) const;
    void save(string name) const;
    void load(string name);
    void save(ofstream& out) const;
    void load(ifstream& in);
    
    bool operator==(const Boundary& bnd) const;
    bool operator!=(const Boundary& bnd) const;
    bool operator%(const Boundary& bnd) const;
    
protected:
    void fourier_interpolation();
    Matrix computeFcoeff(Matrix& Gamma, Real tol, int Nmin);
    Vector FourierInterpolation(Vector curve, Real tol, int Nmin);
    void find_closest_point(const Point& p, Real& t0, Real& err) const;
    void find_closest_discretized_point(const Real t, Real& tj, unsigned& j,const Vector& dtc) const;
    Vector interpolate_nonsmooth(const Vector& v, const Boundary& bnd) const;
    Vector interpolate_smooth(const Vector& v, const Boundary& bnd) const;
    Vector interpolate_smooth_new(const Vector& v, const Boundary& bnd) const;
    
    unsigned n_pts;    
    Real h;
    Real length;
    Vector t;
    Vector dt;
    Matrix G;
    Matrix dG;
    Matrix ddG;
    Matrix N;
    Vector norm_dG;
    Vector norm_ddG;
    Vector measure;
    
    bool smooth_boundary;
    FSeries FGx, FGy;
    FSeries dFGx, dFGy;
    FSeries ddFGx, ddFGy;
    
    vector<array<unsigned,3>> sides_dofs;
    vector<array<unsigned,2>> corners;
    
//    vector<array<unsigned,2>> side_to_dofs_map;
};

}// end namespace uBidomain

ostream& operator<<(ostream& os, const uBidomain::Boundary& bnd);

#endif /* BOUNDARY_H */

