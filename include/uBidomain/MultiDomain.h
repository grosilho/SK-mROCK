#ifndef MULTIDOMAIN_H
#define MULTIDOMAIN_H

#include "Domain.h"
#include "Parameters.h"

typedef vector<array<unsigned,2>> SidesToBoundaryMap;
typedef vector<array<int,6>> GlobalToLocalMap;

namespace uBidomain
{
    
class MultiDomain
{
public:
    MultiDomain();
    MultiDomain(const Boundary& outer_bnd_, //the outer boundary
                const vector<Side>& sides_,//list of internal sides used to construct the inner domains
                const vector<SidesToBoundaryMap> inner_dom_desc);//this tells how to construct the inner domains from the list of internal sides
    
    void get_multidomain_desc(const int n, Parameters& param,
                              Boundary& outer_bnd, vector<Side>& sides, vector<SidesToBoundaryMap>& inner_dom_desc);
    
    void init_MultiDomain(const Boundary& outer_bnd_, const vector<Side>& sides_, 
                          const vector<SidesToBoundaryMap> inner_dom_desc);
    
    void add_AiTViAi(Matrix& F, const Matrix& Vi, const Real alpha, unsigned i, bool sign_change) const;
    void add_AiVAiT(Matrix& Vi, const Matrix& F, const Real alpha, unsigned i, bool sign_change) const;
    void add_AiVAgT(Matrix& Vi, const Matrix& V, const Real alpha, unsigned i, bool sign_change) const;
    void add_AgVAiT(Matrix& Vi, const Matrix& V, const Real alpha, unsigned i, bool sign_change) const;
    void add_AgVAgT(Matrix& Vi, const Matrix& V, const Real alpha) const;
    Matrix global_to_local_left(const Matrix& V, unsigned i, bool sign_change) const;
    Matrix global_to_local_right(const Matrix& V, unsigned i, bool sign_change) const;
    Matrix local_to_global_left(const Matrix& Vi, unsigned i, bool sign_change) const;
    Matrix local_to_global_right(const Matrix& Vi, unsigned i, bool sign_change) const;
    Matrix global_to_gap_junct_left(const Matrix V) const;
    Matrix global_to_gap_junct_right(const Matrix V) const;
    Matrix gap_junct_to_global_left(const Matrix V) const;
    Matrix gap_junct_to_global_right(const Matrix V) const;
    
    Vector remove_measure_gap_junct(const Vector& Vg) const;
    Vector add_measure_gap_junct(const Vector& Vg) const;
    void remove_free_constant(vector<Vector>& v_vec) const;
    
    vector<Matrix> get_PS(vector<bool>& repeated) const;
    
    template<class Sol>
    Vector eval_trace_globally(const Sol& s);
    
    template<class Sol>
    Vector eval_trace_locally(const Sol& s, unsigned i);
    template<class Sol>
    Vector eval_normal_derivative_locally(const Sol& s, unsigned i);
    template<class Sol>
    Vector eval_trace_gap_junct(const Sol& s);
    
    void mesh();
    
    unsigned get_N_domains() const;
    const Domain& get_domain(unsigned i) const;
    Matrix get_coords(unsigned i) const;
    Matrix get_coords_gap_junct() const;
    Vector get_angle_gap_junct() const;
    Vector get_kappa_gap_junct() const;
    unsigned get_n_dofs() const;
    unsigned get_n_outer_bnd_dofs() const;
    unsigned get_n_transmembrane_dofs() const;
    
    ostream& write(ostream& os, string name);
    
    void compute_boundary_errors(const vector<Vector>& traces, const vector<Vector>& derivs,
                                 const vector<Domain>& ref_doms, 
                                 const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                                 Array& err_traces, Array& err_derivs) const;
    void compute_boundary_errors_inner(const vector<Vector>& traces, const vector<Vector>& derivs,
                                 const vector<Domain>& ref_doms, 
                                 const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                                 Array& err_traces, Array& err_derivs) const;
    void compute_boundary_errors(const vector<Vector>& traces, const vector<Vector>& derivs,                                  
                                 const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                                 Array& err_traces, Array& err_derivs) const;
    void compute_boundary_errors_inner(const vector<Vector>& traces, const vector<Vector>& derivs,                                  
                                 const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                                 Array& err_traces, Array& err_derivs) const;
    Real compute_domain_error(const Domain& dom, const Vector& trace, const Vector& deriv,
                              const Domain& ref_dom, const Vector& ref_trace, const Vector& ref_deriv) const;    
    void compute_domain_errors(const vector<Vector>& traces, const vector<Vector>& derivs,
                               const vector<Domain>& ref_doms, 
                               const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                               Array& err_doms) const;
    void compute_domain_errors(const vector<Vector>& traces, const vector<Vector>& derivs,
                               const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                               Array& err_doms) const;
   
   
    void save_domains(const string ref_name) const;
    vector<Domain> load_domains(const string ref_name) const;
    void save_boundary_data(const string ref_name, const vector<Vector>& traces, const vector<Vector>& derivs) const;
    void load_boundary_data(const string ref_name, vector<Vector>& traces, vector<Vector>& derivs) const;
    
    void save(const string ref_name, const vector<Vector>& traces, const vector<Vector>& derivs) const;
    void load(const string ref_name, vector<Vector>& traces, vector<Vector>& derivs, vector<Domain>& domains) const;
    
protected:
    unsigned N_sides;
    Boundary outer_bnd;
    vector<Side> sides;
    unsigned N_domains;
    unsigned N_inner_domains;
    vector<Domain> doms;
    vector<GlobalToLocalMap> doms_GtL_maps;
    vector<array<unsigned,2>> global_dofs;
    GlobalToLocalMap gap_junct_GtL_map;
    vector<unsigned> N_local_dofs;
    unsigned N_global_dofs;
    unsigned N_gap_junct_dofs;
    unsigned N_transmembrane_dofs;
    unsigned N_outer_bnd_dofs;
    Real h_min;
    Real h_max;
        
    Matrix FG;    
    
    Parameters param;
};
    
}

template<class Sol>
Vector uBidomain::MultiDomain::eval_trace_globally(const Sol& s)
{
    Vector trace = Vector::Zero(N_global_dofs);
    
    for(auto side=0u;side<sides.size();side++)
    {
        for(auto i=0u;i<global_dofs[side][1];i++)
            trace[global_dofs[side][0]+i] = s.eval(sides[side].get_G()(i,0),
                                                   sides[side].get_G()(i,1));
    }
    
    return trace;
}

template<class Sol>
Vector uBidomain::MultiDomain::eval_trace_locally(const Sol& s, unsigned i)
{
    return doms[i].eval_trace_on_bnd(s);
}

template<class Sol>
Vector uBidomain::MultiDomain::eval_normal_derivative_locally(const Sol& s, unsigned i)
{
    return doms[i].eval_normal_derivative_on_bnd(s);
}

template<class Sol>
Vector uBidomain::MultiDomain::eval_trace_gap_junct(const Sol& s)
{
    Vector trace(N_gap_junct_dofs);
    
    for(auto j=0u;j<gap_junct_GtL_map.size();j++)
    {
        auto side = gap_junct_GtL_map[j];
        for(auto i=0u;i<side[3];i++)
            trace(side[2]+i) = s.eval(sides[side[0]].get_G()(i,0),sides[side[0]].get_G()(i,1));
    }
    
    return trace;
}

#endif /* MULTIDOMAIN_H */

