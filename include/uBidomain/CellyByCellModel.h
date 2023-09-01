
#ifndef CELLYBYCELLMODEL_H
#define CELLYBYCELLMODEL_H

#include "MultiDomain.h"
#include "Parameters.h"

//typedef Eigen::FullPivHouseholderQR<Matrix> LinearSolver;
//typedef Eigen::LDLT<Matrix> LinearSolver;
typedef Eigen::PartialPivLU<Matrix> LinearSolver;

namespace uBidomain
{
    
class CellByCellModel: public MultiDomain
{
public:
    CellByCellModel(const Boundary& outer_bnd_, 
                const vector<Side>& sides_,
                const vector<SidesToBoundaryMap> inner_dom_desc);
    CellByCellModel(MultiDomain mdom);
    CellByCellModel(const int n, Parameters& param_);
    
    Vector psi(const Vector& V);
    Vector psi(const Vector& V, unsigned i);
    vector<Vector> get_u(const Vector& V) const;
    vector<Vector> get_du(const Vector& V) const;
    void get_u_du(const Vector& V, vector<Vector>& u, vector<Vector>& du) const;
    void add_outer_bnd_data(Vector& trace, Vector& deriv) const;
    
    Real get_sigma(unsigned i) const;
    Real get_Cm() const;
    Matrix get_VtoLambdaM_nomeasure() const;
    
    void get_boundary_data(const Vector& V, vector<Vector>& traces, vector<Vector>& derivs) const;
    void write_sol_and_mesh(const string sol_name, const unsigned n, 
                            const Vector& V) const;
    
    void save_reference_solution(const string ref_name, const Vector& V) const;
    void compute_errors(string ref_name, const Vector& V, Array& err_doms, Array& err_traces, Array& err_derivs) const;
    
    void write_errors(const string filename, const Real dt, const Array& err_doms, const Array& err_traces, const Array& err_derivs) const;
    void write_errors(const string filename, const vector<string> vars_names, const vector<Real> vars) const;
    void print_error(const Array& errors) const;
    
protected:
    void init_parameters();
    void build_system();
    
    
    vector<Real> sigma;
    Real kappa_t;
    Real kappa_l;
    Real Cm;
    
    Matrix F;
    Matrix G;
    Matrix A;
    vector<Vector> ones;
    vector<Matrix> PS;
    vector<Matrix> invPSp;
    
    Vector rhs;
    LinearSolver solver;
    
    bool with_helper_matrices;
    Matrix VtoLambdaM_nomeasure;
    Matrix VtoLambdaM;
    Matrix VtoLambdaG;
    Matrix VtoBeta;
    
    bool artificially_symmetrize;
    bool exact_kernel;
    
    Parameters param;
};

template<typename T>
std::ofstream& operator<<(std::ofstream& out, const vector<T>& v)
{
    for(auto i=0u;i<v.size()-1;i++)
        out<<v[i]<<", ";
    out<<v.back();
    
    return out;
}

}

#endif /* CELLYBYCELLMODEL_H */

