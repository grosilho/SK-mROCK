#include <iomanip>

#include "CellyByCellModel.h"
#include "Utils.h"
#include "Side.h"
#include <Eigen/Eigenvalues>

using namespace uBidomain;

CellByCellModel::CellByCellModel(const Boundary& outer_bnd_, 
                const vector<Side>& sides_,
                const vector<SidesToBoundaryMap> inner_dom_desc)
:MultiDomain(outer_bnd_,sides_,inner_dom_desc)
{
    init_parameters();
}

CellByCellModel::CellByCellModel(MultiDomain mdom)
:MultiDomain(mdom)
{
    init_parameters();
}

CellByCellModel::CellByCellModel(const int n, Parameters& param_)
:param(param_)
{    
    init_parameters();

    Boundary outer_bnd_;
    vector<Side> sides_;
    vector<SidesToBoundaryMap> inner_dom_desc;

    get_multidomain_desc(n,param,outer_bnd_,sides_,inner_dom_desc);
    init_MultiDomain(outer_bnd_,sides_,inner_dom_desc);
        
    sigma.clear();
    sigma.resize(N_domains,param.P20_si);
    sigma[0] = param.P20_se;
    
    build_system();
}

void CellByCellModel::init_parameters()
{
    with_helper_matrices = true;
    artificially_symmetrize = true;
    
    // Values from Silvio Weidman, J. Physiol 210 (1970).
//    Real Rm = 9.1; //kOhm*cm^2
//    kappa = 1./Rm; //mS/cm^2
//    Real Ri = 0.470; //kOhm*cm
//    Real sigma_inner = 1./Ri; //mS/cm
//    Real Re = 0.047; //kOhm*cm
//    Real sigma_outer = 1./Re; //mS/cm
//    Cm = 0.81; //uF/cm^2
    
    // Values from Tveito, Jaeger, Kuchta et al 2017
//    Cm = 1.0;
//    Real sigma_inner = 5.;
//    Real sigma_outer = 20.;
//    Real Rm = 0.0015;
//    kappa = 1./Rm;
    
    //Values from Stinstra, Hopenfeld, MacLeod, 2005
    Cm = 1.0; //uF/cm^2
    Real sigma_inner = 3.; // mS/cm // 3 for longitudinal, 1.5 for trasversal
    Real sigma_outer = 20.;
    Real Rm_l = 0.00145; //kOhm*cm^2
    Real Rm_t = 0.00145; 
    
    // update with command line options
    if(param.P20_si>0.)
        sigma_inner = param.P20_si;
    else
        param.P20_si = sigma_inner;
    if(param.P20_se>0.)
        sigma_outer = param.P20_se;
    else
        param.P20_se = sigma_outer;
    if(param.P20_Rl>0.)
        Rm_l = param.P20_Rl;
    else
        param.P20_Rl = Rm_l;
    if(param.P20_Rt>0.)
        Rm_t = param.P20_Rt;
    else
        param.P20_Rt = Rm_t;
    
    kappa_l = 1./Rm_l;
    kappa_t = 1./Rm_t;        
    
//    ofstream out(string("tmp_geo.m"), ios::out);
//    write(out,"mdom");
//    out.close();      
}

void CellByCellModel::build_system()
{
    F = Matrix::Zero(N_global_dofs,N_global_dofs);
    G = Matrix::Zero(N_global_dofs,N_inner_domains);
    A = Matrix::Zero(N_transmembrane_dofs+N_gap_junct_dofs+N_inner_domains,
                     N_transmembrane_dofs+N_gap_junct_dofs+N_inner_domains);
    ones.clear();
    ones.resize(N_domains,Vector());
    
    cout<<"Building CellByCell model for:"<<endl;
    cout<<"N transmembrane dofs = "<<N_transmembrane_dofs<<endl;
    cout<<"N gap junction dofs = "<<N_gap_junct_dofs<<endl;
    cout<<"N inner domains = "<<N_inner_domains<<endl;
    
    cout<<"Computing PS..."<<flush;
    vector<bool> repeated;
    PS = get_PS(repeated);
    cout<<"done. "<<flush;    
    
    if(artificially_symmetrize)
    {
        cout<<"Symmetrize..."<<flush;
        for(auto& PSi:PS)
            PSi = (PSi+PSi.transpose().eval())/2.;
        cout<<"done. "<<flush;
    }
    
    cout<<"Populate G..."<<flush;
    for(unsigned i=0;i<N_domains;i++)
        ones[i] = Vector::Ones(N_local_dofs[i]);

    for(unsigned i=0;i<N_inner_domains;i++)
        G.col(i) = local_to_global_left(ones[i+1],i+1,true);
    cout<<"done. "<<flush;
    
    cout<<"Populate F..."<<flush;
    invPSp.clear();
    invPSp.resize(N_domains);
    vector<Matrix> PSp(N_domains);
    vector<Real> alpha(N_domains);
    for(unsigned i=0;i<N_domains;i++)
    {    
        if(repeated[i])
        {
            alpha[i]=alpha[i-1];
            PSp[i]=PSp[i-1];
            invPSp[i]=invPSp[i-1];
        }
        else
        {
            cout<<"inverting matrix PSp..."<<flush;
            alpha[i] = 1./pow(PS[i].cols(),1.5);
            PSp[i] = PS[i]+alpha[i]*ones[i]*ones[i].transpose();
            invPSp[i] = PSp[i].inverse();            
        }
    }
    cout<<"done. "<<flush;
    
    cout<<"Projections..."<<flush;
    for(unsigned i=0;i<N_domains;i++)
        add_AiTViAi(F, invPSp[i], -(1./sigma[i]), i, true);        
    cout<<"done. "<<flush;
    
    cout<<"Build matrix blocks F00,F0g,..."<<flush;
    Matrix A11 = Matrix::Zero(N_transmembrane_dofs,N_transmembrane_dofs);
    add_AiVAiT(A11, F, 1., 0, false);
    Matrix A12 = Matrix::Zero(N_transmembrane_dofs,N_gap_junct_dofs);
    add_AiVAgT(A12, F, 1., 0, false);    
    Matrix A21 = Matrix::Zero(N_gap_junct_dofs,N_transmembrane_dofs);
    add_AgVAiT(A21, F, 1., 0, false); 
     
    Vector kappa = get_kappa_gap_junct();
    Matrix A22 = Matrix::Zero(N_gap_junct_dofs,N_gap_junct_dofs);
    A22.diagonal()= -(1./kappa.array());
    add_AgVAgT(A22, F, 1.);
    
    Matrix A13 = global_to_local_left(G,0,false);
    Matrix A23 = global_to_gap_junct_left(G);
    Matrix A31 = local_to_global_right(G.transpose().eval(),0,false);
    Matrix A32 = gap_junct_to_global_right(G.transpose().eval());
    Matrix A33 = Matrix::Zero(N_inner_domains,N_inner_domains);
    cout<<"done. "<<flush;
    
    if(!with_helper_matrices)
    {
        cout<<"Populate A..."<<flush;
        //First line
        A.block(0,0,N_transmembrane_dofs,N_transmembrane_dofs)
                = A11;
        A.block(0,N_transmembrane_dofs,N_transmembrane_dofs,N_gap_junct_dofs)
                = A12;
        A.block(0,N_transmembrane_dofs+N_gap_junct_dofs,N_transmembrane_dofs,N_inner_domains)
                = A13;
        //Second line
        A.block(N_transmembrane_dofs,0,N_gap_junct_dofs,N_transmembrane_dofs)
                = A21;
        A.block(N_transmembrane_dofs,N_transmembrane_dofs,N_gap_junct_dofs,N_gap_junct_dofs)
                = A22;
        A.block(N_transmembrane_dofs,N_transmembrane_dofs+N_gap_junct_dofs,N_gap_junct_dofs,N_inner_domains)
                = A23;
        //Third line
        A.block(N_transmembrane_dofs+N_gap_junct_dofs,0,N_inner_domains,N_transmembrane_dofs)
                = A31;
        A.block(N_transmembrane_dofs+N_gap_junct_dofs,N_transmembrane_dofs,N_inner_domains,N_gap_junct_dofs)
                = A32;
        A.block(N_transmembrane_dofs+N_gap_junct_dofs,N_transmembrane_dofs+N_gap_junct_dofs,N_inner_domains,N_inner_domains)
                = A33;    
        cout<<"done. "<<flush;
        
        cout<<"Solve A..."<<flush;
        solver.compute(A);
        rhs = Vector::Zero(N_transmembrane_dofs+N_gap_junct_dofs+N_inner_domains);
        cout<<"done. "<<flush;    
    } 
    else
    {
        Matrix I11, A11inv, K2, K3, K1K2, K1K3, tmp1, tmp2;
        LinearSolver K1;
//        cout<<"Computing helper matrices..."<<flush;
//
//        I11 = Matrix::Identity(A11.rows(),A11.cols());
//        A11inv = A11.partialPivLu().solve(I11);
//
//        K3 = -A21*A11inv;
//        K1.compute( A22-A21*A11inv*A12 );
//        K2 = A21*(A11inv*A13)-A23;        
//        K1K2 = K1.solve(K2);
//        K1K3 = K1.solve(K3);
//
//        tmp1 = (A32-A31*A11inv*A12)*K1K2-A31*A11inv*A13;
//        tmp2 = A31*A11inv*(A12*K1K3-I11)-A32*K1K3;
//        VtoBeta = tmp1.partialPivLu().solve(tmp2);
//
//        VtoLambdaM = A11inv*((I11-A12*K1K3)-(A13+A12*K1K2)*VtoBeta);
//        VtoLambdaG = K1K3+K1K2*VtoBeta;        
//        VtoLambdaM_nomeasure = Matrix::Zero(A11.rows(),A11.cols());
//        for(unsigned i=0;i<A11.cols();i++)
//            VtoLambdaM_nomeasure.col(i) = doms[0].remove_inner_measure(VtoLambdaM.col(i));    
//        cout<<"done."<<flush<<endl;
        
        //alternative way
        cout<<"Computing helper matrices..."<<flush;
        I11 = Matrix::Identity(A11.rows(),A11.cols());
        A11inv = A11.partialPivLu().solve(I11);

        Matrix A11invA13 = A11inv*A13;
        K3 = -A21*A11inv;
        Matrix tmp4 = A22+K3*A12;
        K1.compute( tmp4 );
        K2 = A21*(A11invA13)-A23;        
        K1K2 = K1.solve(K2);        
        K1K3 = K1.solve(K3);
        
        Matrix tmp3 = A11inv*(A12*K1K3)-A11inv; //A11inv*(A12*K1K3-I11)
        tmp1 = (A32+A31*K3.transpose())*K1K2-A31*A11invA13;
        tmp2 = A31*tmp3-A32*K1K3;
        VtoBeta = tmp1.partialPivLu().solve(tmp2);

        VtoLambdaM = -tmp3+(K3.transpose()*K1K2-A11invA13)*VtoBeta;
        VtoLambdaG = K1K3+K1K2*VtoBeta;        
        VtoLambdaM_nomeasure = Matrix::Zero(A11.rows(),A11.cols());
        for(unsigned i=0;i<A11.cols();i++)
            VtoLambdaM_nomeasure.col(i) = doms[0].remove_inner_measure(VtoLambdaM.col(i));    
        cout<<"done."<<flush<<endl;                
        
    // For coarse problems the matrix VtoLambdaM_nomeasure might
    // contain positive eigenvalues, which lead to instability
    // Activate the next portion of code to filter these eigenvalues out
//    cout<<"Computing eigendecomposition"<<endl;
//    Eigen::EigenSolver<Matrix> eigs(VtoLambdaM_nomeasure);
//    cout<<"done"<<endl;
//    Eigen::MatrixXcd eig_vects = eigs.eigenvectors();
//    Eigen::VectorXcd eig_vals = eigs.eigenvalues();
//    unsigned found_pos = 0;
//    for(unsigned i=0;i<eig_vals.size();i++)
//        if(eig_vals[i].real()>0)
//        {
//            found_pos++;
//            eig_vals[i]*=-1.;
//        }
//    if(found_pos>0)
//    {
//        cout<<"Found "<<found_pos<<" positive eigenvalues"<<endl;
//        Eigen::MatrixXcd D(VtoLambdaM_nomeasure.rows(),VtoLambdaM_nomeasure.cols());
//        D.diagonal() = eig_vals;
//        cout<<"computing new matrix"<<endl;
//        Eigen::MatrixXcd newV = eig_vects*D*eig_vects.inverse();
//        cout<<"error = "<<(VtoLambdaM_nomeasure -newV).norm()<<endl;
//        VtoLambdaM_nomeasure = newV.real();
//    }
    
//    ofstream ofile("VtoL.m",ios::out);
//    write_matrix(ofile,VtoLambdaM_nomeasure,"A");
//    ofile.close();
    }
}

Vector CellByCellModel::psi(const Vector& V)
{
    if(!with_helper_matrices)
    {
        rhs.head(N_transmembrane_dofs) = V;
        Vector sol = solver.solve(rhs);    
        return doms[0].remove_inner_measure(sol.head(N_transmembrane_dofs));
    }
    else
        return VtoLambdaM_nomeasure*V;
}

Vector CellByCellModel::psi(const Vector& V, unsigned i)
{
    rhs.head(N_transmembrane_dofs+N_gap_junct_dofs) = V;
    
    Matrix A(F.rows()+G.cols(),F.cols()+G.cols());
    A.block(0,0,F.rows(),F.cols())= F;
    A.block(0,F.cols(),G.rows(),G.cols()) = G;
    A.block(F.rows(),0,G.cols(),G.rows()) = G.transpose();
    A.block(F.rows(),G.rows(),G.cols(),G.cols()) = Matrix::Zero(G.cols(),G.cols());
    
    LinearSolver lin_solv;
    lin_solv.compute(A);
    
    Vector sol = lin_solv.solve(rhs);
    Vector lambda = sol.head(N_transmembrane_dofs+N_gap_junct_dofs);
    
    return -doms[i].remove_inner_measure(global_to_local_left(lambda,i,true));
}

vector<Vector> CellByCellModel::get_u(const Vector& V) const
{
    Vector lambda;
    Vector beta(N_domains);
    beta(0) = 0.;
        
    if(!with_helper_matrices)
    {
        Vector rhs = Vector::Zero(N_transmembrane_dofs+N_gap_junct_dofs+N_inner_domains);
        rhs.head(N_transmembrane_dofs) = V;
        Vector sol = solver.solve(rhs);
        
        Vector lambda_m = sol.head(N_transmembrane_dofs);
        Vector lambda_g = sol.segment(N_transmembrane_dofs,N_gap_junct_dofs);
        lambda = local_to_global_left(lambda_m,0,false)
                        +gap_junct_to_global_left(lambda_g);

        beta.tail(N_inner_domains)= sol.tail(N_inner_domains);
    }
    else
    {
        Vector lambda_m = VtoLambdaM*V;
        Vector lambda_g = VtoLambdaG*V;
        lambda = local_to_global_left(lambda_m,0,false)
                        +gap_junct_to_global_left(lambda_g);
        
        beta.tail(N_inner_domains) = VtoBeta*V;
    }
    
    vector<Vector> u(N_domains);
    for(unsigned i=0;i<N_domains;i++)
    {
        u[i] = -(1./sigma[i])*invPSp[i]*global_to_local_left(lambda,i,true)
                +beta[i]*ones[i];
    }
    
    return u;
}

vector<Vector> CellByCellModel::get_du(const Vector& V) const
{
    Vector lambda;
    if(!with_helper_matrices)
    {
        Vector rhs = Vector::Zero(N_transmembrane_dofs+N_gap_junct_dofs+N_inner_domains);
        rhs.head(N_transmembrane_dofs) = V;
        Vector sol = solver.solve(rhs);
        
        Vector lambda_m = sol.head(N_transmembrane_dofs);
        Vector lambda_g = sol.segment(N_transmembrane_dofs,N_gap_junct_dofs);
        lambda = local_to_global_left(lambda_m,0,false)
                        +gap_junct_to_global_left(lambda_g);
    }
    else
    {
        Vector lambda_m = VtoLambdaM*V;
        Vector lambda_g = VtoLambdaG*V;
        lambda = local_to_global_left(lambda_m,0,false)
                        +gap_junct_to_global_left(lambda_g);
    }
    
    vector<Vector> du(N_domains);
    for(unsigned i=0;i<N_domains;i++)
    {        
        du[i] = -(1./sigma[i])*doms[i].remove_inner_measure(global_to_local_left(lambda,i,true));
    }

    return du;
}

void CellByCellModel::get_u_du(const Vector& V, vector<Vector>& u, vector<Vector>& du) const
{
    Vector lambda;
    Vector beta(N_domains);
    beta(0) = 0.;
        
    if(!with_helper_matrices)
    {
        Vector rhs = Vector::Zero(N_transmembrane_dofs+N_gap_junct_dofs+N_inner_domains);
        rhs.head(N_transmembrane_dofs) = V;
        Vector sol = solver.solve(rhs);
        
        Vector lambda_m = sol.head(N_transmembrane_dofs);
        Vector lambda_g = sol.segment(N_transmembrane_dofs,N_gap_junct_dofs);
        lambda = local_to_global_left(lambda_m,0,false)
                        +gap_junct_to_global_left(lambda_g);

        beta.tail(N_inner_domains)= sol.tail(N_inner_domains);
    }
    else
    {
        Vector lambda_m = VtoLambdaM*V;
        Vector lambda_g = VtoLambdaG*V;
        lambda = local_to_global_left(lambda_m,0,false)
                        +gap_junct_to_global_left(lambda_g);
        
        beta.tail(N_inner_domains) = VtoBeta*V;
    }
    
    u.resize(N_domains);
    du.resize(N_domains);
    for(unsigned i=0;i<N_domains;i++)
    {
        u[i] = -(1./sigma[i])*invPSp[i]*global_to_local_left(lambda,i,true)
                +beta[i]*ones[i];
        du[i] = -(1./sigma[i])*doms[i].remove_inner_measure(global_to_local_left(lambda,i,true));
    }

}

void CellByCellModel::add_outer_bnd_data(Vector& trace, Vector& deriv) const
{
    Vector trace_ext(N_outer_bnd_dofs+N_local_dofs[0]);
    Vector deriv_ext = Vector::Zero(N_outer_bnd_dofs+N_local_dofs[0]);
    trace_ext.tail(N_local_dofs[0]) = trace;
    deriv_ext.tail(N_local_dofs[0]) = deriv;
    trace_ext.head(N_outer_bnd_dofs) = doms[0].compute_outer_trace(trace,deriv_ext.head(N_outer_bnd_dofs));
    
    trace = std::move(trace_ext);
    deriv = std::move(deriv_ext);
}

Real CellByCellModel::get_sigma(unsigned i) const
{
    return sigma[i];
}

Real CellByCellModel::get_Cm() const
{
    return Cm;
}

Matrix CellByCellModel::get_VtoLambdaM_nomeasure() const
{
    if(!with_helper_matrices)
        cout<<"VtoLambdaM_nomeasure is probably not initialized"<<endl;
    
    return VtoLambdaM_nomeasure;
}

void CellByCellModel::get_boundary_data(const Vector& V,
                                vector<Vector>& traces, vector<Vector>& derivs) const
{
    get_u_du(V,traces,derivs);
    add_outer_bnd_data(traces[0],derivs[0]);
}

void CellByCellModel::write_sol_and_mesh(const string sol_name, const unsigned n, 
                                         const Vector& V) const
{
    string dom = "_dom_";
    
    vector<Vector> traces;
    vector<Vector> derivs;
    get_boundary_data(V,traces,derivs);
    
    for(unsigned i=0;i<N_domains;i++)
    {
        vector<Real> ui = doms[i].eval_on_mesh(traces[i],derivs[i]);
        doms[i].write_sol_and_mesh(sol_name+dom+to_string(i)+"_n_"+to_string(n),ui);
    }
}

void CellByCellModel::save_reference_solution(const string ref_name, const Vector& V) const
{    
    vector<Vector> traces;
    vector<Vector> derivs;
    get_boundary_data(V,traces,derivs);
    
    save(ref_name,traces,derivs);
}

vector<Domain> MultiDomain::load_domains(const string ref_name) const
{
    vector<Domain> ref_doms;
    ifstream ifile(ref_name+"_domains.bin", std::ios::in);
    for(unsigned i=0;i<N_domains;i++)
        ref_doms.push_back(Domain(ifile));
    ifile.close();
    return ref_doms;
}

void MultiDomain::load_boundary_data(const string ref_name, vector<Vector>& traces, vector<Vector>& derivs) const
{
    ifstream ifile(ref_name+"_solutions.bin", std::ios::in);
    traces.resize(N_domains);
    derivs.resize(N_domains);
    
    for(unsigned i=0;i<N_domains;i++)
    {
        read_binary(ifile, traces[i]);
        read_binary(ifile, derivs[i]);
    }
    ifile.close();
}

void CellByCellModel::compute_errors(string ref_name, const Vector& V, 
        Array& err_doms, Array& err_traces, Array& err_derivs) const
{
    vector<Domain> ref_doms;    
    vector<Vector> ref_traces;
    vector<Vector> ref_derivs;
    load(ref_name,ref_traces,ref_derivs,ref_doms);
    
    vector<Vector> traces;
    vector<Vector> derivs;
    get_boundary_data(V,traces,derivs);      
       
    compute_domain_errors(traces,derivs,
                          ref_doms,ref_traces,ref_derivs,
                          err_doms);   

    compute_boundary_errors(traces,derivs,
            ref_doms,ref_traces,ref_derivs,err_traces,err_derivs);
}

void CellByCellModel::write_errors(const string filename, const Real dt, const Array& err_doms, const Array& err_traces, const Array& err_derivs) const
{
    vector<string> vars_names;
    vector<Real> vars;
    
    vars_names.push_back("dt");
    vars.push_back(dt);
    
    vars_names.push_back("h");
    vars.push_back(h_max);
    
    vars_names.push_back(string("max_dom_err"));
    vars.push_back(err_doms.maxCoeff());

    vars_names.push_back(string("max_trace_err"));
    vars.push_back(err_traces.maxCoeff());

    vars_names.push_back(string("max_deriv_err"));
    vars.push_back(err_derivs.maxCoeff());
    
    for(unsigned i=0;i<N_domains;i++)
    {
        vars_names.push_back(string("dom_err_")+to_string(i));
        vars.push_back(err_doms(i));
        
        vars_names.push_back(string("trace_err_")+to_string(i));
        vars.push_back(err_traces(i));
        
        vars_names.push_back(string("deriv_err_")+to_string(i));
        vars.push_back(err_derivs(i));
    }
    
    write_errors(filename,vars_names,vars);
}

void CellByCellModel::write_errors(const string filename, const vector<string> vars_names, const vector<Real> vars) const
{
    ifstream ifile(filename+".csv");
    bool exists = ifile.good();
    ifile.close();

    ofstream ofile(filename+".csv", std::ios::out | std::ios::app);
    if(!exists)
        ofile<<vars_names<<endl;

    ofile<<std::fixed<<std::setprecision(16);
    ofile<<vars<<endl;
    ofile.close();
}

void CellByCellModel::print_error(const Array& errors) const
{
    // printing results 
    unsigned int prec = 4;
    unsigned int num_size = prec+6;
    unsigned int cell_size = num_size+1;
    unsigned int N = errors.size();
    unsigned int table_width = 1+(cell_size+1)*N;
        
    cout<<setprecision(prec)<<scientific;
    
    string text = "Errors in l2 norm";
    unsigned len = text.size();
    unsigned a= (table_width-len)/2;
    unsigned b= (table_width-len)%2;
    
    cout<<setfill('-')<<setw(a)<<"-"<<text<<setfill('-')<<setw(a+b)<<"-"<<endl;

    cout<<"|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<("Dom "+to_string(i))<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    cout<<"|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<errors(i)<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
}
