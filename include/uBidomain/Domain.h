#ifndef DOMAIN_H
#define DOMAIN_H

#include <memory>
#include <vector>
#include <map>
#include <set>

#include "Boundary.h"

namespace uBidomain
{
class Domain
{
public:
    Domain(const Boundary& bnd_);
    Domain(const vector<Boundary>& bnds_);
    
    Domain(string filename);
    Domain(ifstream& in);

    Matrix SingleLayerOp() const;
    Matrix SingleLayerOp(unsigned bnd1, unsigned bnd2) const;
    Matrix SingleLayerOp(unsigned bnd) const;
    Matrix DoubleLayerOp() const;
    Matrix DoubleLayerOp(unsigned bnd1, unsigned bnd2) const;

    Matrix PoincareSteklow() const;
    void MyPoincareSteklow(Matrix& A, Matrix& B) const;
    Matrix MyPoincareSteklow() const;
    Vector compute_inner_deriv(const Vector& tracei, const Vector& derivo) const;
    Vector compute_outer_trace(const Vector& tracei, const Vector& derivo) const;
    void compute_helper_matrices();
    Vector DirichletToNeumann(const Vector& trace) const;
    Vector NeumannToDirichlet(const Vector& deriv) const;
    
    void mesh();
    void write_sol_and_mesh(const string sol_name, const vector<Real>& u) const;

    template<class Sol>
    vector<Real> eval_on_mesh(const Sol& s) const;

    Real eval_on_point(const Point p, const Vector& trace, const Vector& deriv) const;

    vector<Real> eval_on_points(const vector<Real>& pts, const Vector& trace, const Vector& deriv) const;
    
    vector<Real> eval_on_mesh(const Vector& trace, const Vector& deriv) const;

    vector<Real> get_non_boundary_points() const;
    
    template<class Sol>
    Vector eval_trace_on_bnd(const Sol& s) const;
    template<class Sol>
    Vector eval_normal_derivative_on_bnd(const Sol& s) const;

    Vector remove_measure(const Vector& deriv) const;
    Vector add_measure(const Vector& deriv) const;
    Vector remove_measure(const Vector& deriv, unsigned bnd) const;
    Vector add_measure(const Vector& deriv, unsigned bnd) const;
    Vector remove_inner_measure(const Vector& deriv) const;
    Vector add_inner_measure(const Vector& deriv) const;
    Vector remove_outer_measure(const Vector& deriv) const;
    Vector add_outer_measure(const Vector& deriv) const;
    
    Vector get_outer(const Vector& v) const;
    Vector get_inner(const Vector& v) const;        

    Real L2norm(const Vector& v) const;
    Real L2norm_inner(const Vector& v) const;
    Real L2norm(const Vector& v, unsigned i) const;
    Real Linfnorm(const Vector& v) const;
    
    Vector interpolate(const Vector& v, const Domain& dom) const;
    Vector interpolate_inner(const Vector& v, const Domain& dom) const;
    
    Point point(unsigned ind) const;
    const vector<Boundary>& get_bnds() const;
    unsigned get_tot_n_dofs() const;
    Real get_h_min() const;
    Real get_h_max() const;
    
    ostream& write_boundary_sol_and_coords(ostream& os, const Vector& v, string name) const;
    
    ostream& write(ostream& os, string name) const;
    void save(string name) const;
    void load(string name);
    void save(ofstream& out) const;
    void load(ifstream& in);
    
    bool operator==(const Domain& dom) const;
    bool operator%(const Domain& dom) const;

protected:
    void initialize();
    void compute_h_min_max();
    bool is_point_on_segment(const Point& p, const Point& q, 
                             const Point& x, Real& t, Real& s);
    
    unsigned N_bnds;
    unsigned N_inner_bnds;
    unsigned N_outer_dofs;
    unsigned N_inner_dofs;
    Real h_min;
    Real h_max;
    
    vector<Boundary> bnds;
    
    vector<unsigned> bnd_dofs;
    vector<unsigned> bnd_size;

    bool meshed;
    vector<Real> points;
    vector<int> elems;
    vector<int> is_point_on_bnd;
    vector<tuple<unsigned,unsigned,unsigned,Real,Real>> map_to_bnd_points;
    vector<int> boundary_pts;
    vector<int> non_boundary_pts;
    
    Eigen::PartialPivLU<Matrix> V_II_fact;
//    Eigen::LDLT<Matrix> V_II_fact;
    
//    Eigen::ColPivHouseholderQR<Matrix> C_fact;
    Eigen::PartialPivLU<Matrix> C_fact;
    
    Matrix V_OO,V_OI,V_IO,V_II,K_OO,K_OI,K_IO,K_II,I_O,I_I;
    Matrix matD;
};

}// end namespace uBidomain

template<class Sol>
vector<Real> uBidomain::Domain::eval_on_mesh(const Sol& s) const
{
    if(points.empty() || elems.empty())
    {
        cerr<<"ERROR: Run Domain::mesh() first."<<endl;
        return vector<Real>();
    }
    
    vector<Real> u(points.size()/2,0.);
    for(auto i=0u;i<u.size();i++)
        u[i] = s.eval(points[2*i],points[2*i+1]);
    return u;
}

template<class Sol>
Vector uBidomain::Domain::eval_trace_on_bnd(const Sol& s) const
{
    Vector trace(N_outer_dofs+N_inner_dofs);
    
    for(unsigned bnd=0;bnd<N_bnds;bnd++)
    {
        const Matrix& G = bnds[bnd].get_G();
        for(unsigned i=0;i<bnds[bnd].get_n_pts();i++)
            trace(bnd_dofs[bnd] + i) = s.eval(G(i,0),G(i,1));
    }    
    return trace;
}

template<class Sol>
Vector uBidomain::Domain::eval_normal_derivative_on_bnd(const Sol& s) const
{
    Vector deriv(N_outer_dofs+N_inner_dofs);
    
    for(unsigned bnd=0;bnd<N_bnds;bnd++)
    {
        const Matrix& G = bnds[bnd].get_G();
        unsigned n_pts = bnds[bnd].get_n_pts();
        for(unsigned i=0;i<n_pts;i++)
            deriv(bnd_dofs[bnd] + i) = bnds[bnd].get_N(i)*s.eval_grad(G(i,0),G(i,1));        
    }    
    return deriv;
}

ostream& operator<<(ostream& os, const uBidomain::Domain& dom);


#endif /* DOMAIN_H */
