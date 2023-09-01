#include <limits>

#include "gmsh.h"

#include "Domain.h"
#include "lean_vtk.h"
#include "getV1.h"
#include "Utils.h"


using namespace uBidomain;

Domain::Domain(const Boundary& bnd_)
:bnds{bnd_}
{    
    initialize();
    compute_helper_matrices();
}

Domain::Domain(const vector<Boundary>& bnds_)
:bnds{bnds_}
{   
    initialize();
    compute_helper_matrices();
}

Domain::Domain(string filename)
{
    load(filename);
    compute_helper_matrices();
}

Domain::Domain(ifstream& in)
{
    load(in);
    compute_helper_matrices();
}

void Domain::initialize()
{
    N_bnds = bnds.size();
    N_inner_bnds = N_bnds-1;
        
    N_outer_dofs = bnds[0].get_n_pts();
    
    bnd_dofs.resize(N_bnds);
    bnd_size.resize(N_bnds);
    bnd_dofs[0]=0;
    bnd_size[0]=bnds[0].get_n_pts();
    for(auto b=1u;b<N_bnds;b++)
    {
        bnd_size[b] = bnds[b].get_n_pts();
        bnd_dofs[b] = bnd_dofs[b-1]+bnd_size[b-1];
    }
    N_inner_dofs = bnd_dofs.back()+bnd_size.back()-N_outer_dofs;
    
    meshed=false;
    
    compute_h_min_max();
}

void Domain::compute_h_min_max()
{
    h_min = std::numeric_limits<Real>::max();
    h_max = 0.;
    
    for(auto bnd:bnds)
    {
        h_min = min(h_min,bnd.get_h());
        h_max = max(h_max,bnd.get_h());
    }
    
}

Matrix Domain::SingleLayerOp() const
{
    Matrix V(N_outer_dofs+N_inner_dofs,N_outer_dofs+N_inner_dofs);
    
    for(auto i=0u;i<N_bnds;i++)
        for(auto j=0u;j<N_bnds;j++)
            V.block(bnd_dofs[i],bnd_dofs[j],
                    bnd_size[i],bnd_size[j]) = SingleLayerOp(i,j);
    
    return V;
}

Matrix Domain::SingleLayerOp(unsigned bnd1, unsigned bnd2) const
{   
    if(bnd1==bnd2)
        return SingleLayerOp(bnd1);
    
    //if we arrive here the two boundaries are different and we do not care
    //about singularities
    Matrix V(bnd_size[bnd1],bnd_size[bnd2]);
    Eigen::ArrayXd diffx,diffy,normsq;
    
    for(auto j=0u;j<bnd_size[bnd2];j++)
    {
        diffx = bnds[bnd1].get_G().col(0).array()-bnds[bnd2].get_G()(j,0);
        diffy = bnds[bnd1].get_G().col(1).array()-bnds[bnd2].get_G()(j,1);
        normsq = diffx.square()+diffy.square();
        V.col(j) = (-1./4./M_PI)*normsq.log();
    }
    
    return V;
}

Matrix Domain::SingleLayerOp(unsigned bnd) const
{        
    Matrix V(bnd_size[bnd],bnd_size[bnd]);
    Eigen::ArrayXd diffx,diffy,difft,normsq,sinsq;
    
    const Matrix& G = bnds[bnd].get_G();
    const Vector& norm_dG = bnds[bnd].get_norm_dG();
    const Vector& measure = bnds[bnd].get_measure();
    const Vector& t = bnds[bnd].get_t();
    const Vector& dt = bnds[bnd].get_dt();
    
    for(auto j=0u;j<bnd_size[bnd];j++)
    {
        diffx = G.col(0).array()-G(j,0);
        diffy = G.col(1).array()-G(j,1);
        normsq = diffx.square()+diffy.square();
        difft = t.array()-t(j);
        sinsq = (M_PI*difft).sin().square();
        V.col(j) = (normsq/(4.*sinsq)).log();
    }
    V.diagonal() = 2.*(norm_dG.array()/2./M_PI).log();
    
    Matrix V1;
    V1.resize(bnd_size[bnd], bnd_size[bnd]);
    getV1(V1.data(), bnd_size[bnd]);
    
    V = V*(-1./4./M_PI);
    
    for(auto j=0u;j<bnd_size[bnd];j++)
        V1.col(j) = V1.col(j)/dt(j);
    
    V = V+V1;
    
    return V;
}

Matrix Domain::DoubleLayerOp() const
{
    Matrix K(N_outer_dofs+N_inner_dofs,N_outer_dofs+N_inner_dofs);
    
    for(auto i=0u;i<N_bnds;i++)
        for(auto j=0u;j<N_bnds;j++)
            K.block(bnd_dofs[i],bnd_dofs[j],
                    bnd_size[i],bnd_size[j]) = DoubleLayerOp(i,j);
    
    // If a boundary has corners then at the nodes close to a corner the integral is
    // computed incorrectly. We correct it supposing that away from the nodes
    // the integral is more accurate and that (K+I/2)*ones=0 (hence the sum over
    // a line of K must be -0.5).
    for(auto bnd=0u;bnd<N_bnds;bnd++)
    {
        vector<array<unsigned,2>> corners = bnds[bnd].get_corners();
        for(auto corner:corners)
            K(bnd_dofs[bnd]+corner[0],bnd_dofs[bnd]+corner[1])=-0.5-(K.row(bnd_dofs[bnd]+corner[0]).array().sum()-K(bnd_dofs[bnd]+corner[0],bnd_dofs[bnd]+corner[1]));
    }
    
    return K;
}

Matrix Domain::DoubleLayerOp(unsigned bnd1, unsigned bnd2) const
{   
    Matrix K(bnd_size[bnd1],bnd_size[bnd2]);
    Eigen::ArrayXd diffx,diffy,normsq;
    
    const Matrix& G1 = bnds[bnd1].get_G();
    const Matrix& G2 = bnds[bnd2].get_G();
    const Matrix& N2 = bnds[bnd2].get_N();
    const Vector& norm_dG2 = bnds[bnd2].get_norm_dG();
    const Vector& measure2 = bnds[bnd2].get_measure();
    const Vector& dt2 = bnds[bnd2].get_dt();
    
    for(auto j=0u;j<bnd_size[bnd2];j++)
    {
        diffx = G1.col(0).array()-G2(j,0);
        diffy = G1.col(1).array()-G2(j,1);
        normsq = diffx.square()+diffy.square();
        K.col(j) = diffx*N2(j,0)+diffy*N2(j,1);
        K.col(j) = (K.col(j).array()/normsq)/(2.*M_PI);
        K.col(j) = K.col(j)*measure2(j); // comment to remove measure
        
    }
    
    if(bnd1==bnd2)
    {
        const Matrix& ddG2 = bnds[bnd2].get_ddG();
        K.diagonal() = ddG2.col(0).cwiseProduct(N2.col(0))+ddG2.col(1).cwiseProduct(N2.col(1));
        K.diagonal() = (1./2./M_PI)*K.diagonal().cwiseQuotient(2.*norm_dG2).cwiseProduct(dt2); 
        
        //uncomment to remove measure
//        K.diagonal() = K.diagonal().cwiseQuotient(norm_dG2*dt2); 
        
//        for(auto corner:bnds[bnd1].get_corners())
//        {
//            int i = corner[0];
//            int j = corner[1];
//            Real diffx = G1(i,0)-G2(j,0);
//            Real diffy = G1(i,1)-G2(j,1);
//            Real normsq = pow(diffx,2)+pow(diffy,2);
//            Point N = (N2.row(i)+N2.row(j))/2.;
//            K(i,j) = diffx*N(0)+diffy*N(1);
//            K(i,j) = (K(i,j)/normsq)/(2.*M_PI);
//            K(i,j) = K(i,j)*norm_dG2(j)*dt2; // comment to remove measure
//        }
    }
    
    
    return K;
}

Matrix Domain::PoincareSteklow() const
{
    Matrix V = SingleLayerOp();
    Matrix K = DoubleLayerOp();
    Matrix I = Matrix::Identity(V.rows(),V.cols());

    return V.fullPivHouseholderQr().solve(K+0.5*I);
}

void Domain::MyPoincareSteklow(Matrix& A, Matrix& B) const
{
    // If the domain has one boundary (Outer) then A is the Dirichlet to Neumann
    // map and B is empty. If the domain has outer and inner bnd then A,B are such that
    // Av+Bw is the Neumann data on the inner boundary, where v is the dirichlet
    // data on inner bnd and w is Neumann data on outer bnd. 
    
    if(N_inner_bnds>0)
    {
        Matrix V = SingleLayerOp();
        Matrix K = DoubleLayerOp();
        
        Matrix V_OO = V.block(0,0,N_outer_dofs,N_outer_dofs);
        Matrix V_OI = V.block(0,N_outer_dofs,N_outer_dofs,N_inner_dofs);
        Matrix V_IO = V.block(N_outer_dofs,0,N_inner_dofs,N_outer_dofs);
        Matrix V_II = V.block(N_outer_dofs,N_outer_dofs,N_inner_dofs,N_inner_dofs);
        Matrix K_OO = K.block(0,0,N_outer_dofs,N_outer_dofs);
        Matrix K_OI = K.block(0,N_outer_dofs,N_outer_dofs,N_inner_dofs);
        Matrix K_IO = K.block(N_outer_dofs,0,N_inner_dofs,N_outer_dofs);
        Matrix K_II = K.block(N_outer_dofs,N_outer_dofs,N_inner_dofs,N_inner_dofs);
        Matrix I_O = Matrix::Identity(N_outer_dofs,N_outer_dofs);
        Matrix I_I = Matrix::Identity(N_inner_dofs,N_inner_dofs);
        
        Eigen::PartialPivLU<Matrix> V_II_fact;
        // Eigen::ColPivHouseholderQR<Matrix> V_II_fact;
        V_II_fact.compute(V_II);
        
        Eigen::PartialPivLU<Matrix> C_fact;
        // Eigen::ColPivHouseholderQR<Matrix> C_fact;
        C_fact.compute(K_OO+0.5*I_O-V_OI*V_II_fact.solve(K_IO));
        
        Matrix D = V_OI*V_II_fact.solve(K_II+0.5*I_I)-K_OI;
        Matrix E = V_OO-V_OI*V_II_fact.solve(V_IO);
        
        A = V_II_fact.solve(K_II+0.5*I_I+K_IO*C_fact.solve(D));
        B = V_II_fact.solve(K_IO*C_fact.solve(E)-V_IO);
    }
    else
    {
        A = PoincareSteklow();
        B = Matrix();
    }
}

Matrix Domain::MyPoincareSteklow() const
{
    // If the domain has one boundary (Outer) then A is the Dirichlet to Neumann
    // map and B is empty. If the domain has outer and inner bnd then A,B are such that
    // Av+Bw is the Neumann data on the inner boundary, where v is the dirichlet
    // data on inner bnd and w is Neumann data on outer bnd. 
    
    if(N_inner_bnds>0)
    {
//        Matrix V = SingleLayerOp();
//        Matrix K = DoubleLayerOp();
//        
//        Matrix V_OO = V.block(0,0,N_outer_dofs,N_outer_dofs);
//        Matrix V_OI = V.block(0,N_outer_dofs,N_outer_dofs,N_inner_dofs);
//        Matrix V_IO = V.block(N_outer_dofs,0,N_inner_dofs,N_outer_dofs);
//        Matrix V_II = V.block(N_outer_dofs,N_outer_dofs,N_inner_dofs,N_inner_dofs);
//        Matrix K_OO = K.block(0,0,N_outer_dofs,N_outer_dofs);
//        Matrix K_OI = K.block(0,N_outer_dofs,N_outer_dofs,N_inner_dofs);
//        Matrix K_IO = K.block(N_outer_dofs,0,N_inner_dofs,N_outer_dofs);
//        Matrix K_II = K.block(N_outer_dofs,N_outer_dofs,N_inner_dofs,N_inner_dofs);
//        Matrix I_O = Matrix::Identity(N_outer_dofs,N_outer_dofs);
//        Matrix I_I = Matrix::Identity(N_inner_dofs,N_inner_dofs);
        
//        Eigen::PartialPivLU<Matrix> V_II_fact;
//        // Eigen::ColPivHouseholderQR<Matrix> V_II_fact;
//        V_II_fact.compute(V_II);
//        
//        Eigen::PartialPivLU<Matrix> C_fact;
//        // Eigen::ColPivHouseholderQR<Matrix> C_fact;
//        C_fact.compute(K_OO+0.5*I_O-V_OI*V_II_fact.solve(K_IO));
        
        
//        Matrix E = V_OO-V_OI*V_II_fact.solve(V_IO);
        
        return V_II_fact.solve(K_II+0.5*I_I+K_IO*C_fact.solve(matD));
    }
    else
        return PoincareSteklow();
}

//void Domain::V_II_fact_compute(const Matrix& V_II)
//{
//    V_II_fact.compute(V_II);
//}
//
//void Domain::C_fact_compute(const Matrix& tmp)
//{
//    C_fact.compute(tmp);
//}
//
//void Domain::compute_A(const Matrix& V_OO,const Matrix& V_OI,const Matrix& V_IO,const Matrix& V_II,
//                       const Matrix& K_OO,const Matrix& K_OI,const Matrix& K_IO,const Matrix& K_II,
//                       const Matrix& I_O,const Matrix& I_I)
//{
//    matA = V_II_fact.solve(K_II+0.5*I_I+K_IO*C_fact.solve(matD));
//}
//
//void Domain::compute_B(const Matrix& V_OO,const Matrix& V_OI,const Matrix& V_IO,const Matrix& V_II,
//                       const Matrix& K_OO,const Matrix& K_OI,const Matrix& K_IO,const Matrix& K_II,
//                       const Matrix& I_O,const Matrix& I_I)
//{
//    matB = V_II_fact.solve(K_IO*C_fact.solve(matE)-V_IO);
//}
//
//void Domain::compute_D(const Matrix& V_OO,const Matrix& V_OI,const Matrix& V_IO,const Matrix& V_II,
//                       const Matrix& K_OO,const Matrix& K_OI,const Matrix& K_IO,const Matrix& K_II,
//                       const Matrix& I_O,const Matrix& I_I)
//{
//    matD = V_OI*V_II_fact.solve(K_II+0.5*I_I)-K_OI;
//}
//
//void Domain::compute_E(const Matrix& V_OO,const Matrix& V_OI,const Matrix& V_IO,const Matrix& V_II,
//                       const Matrix& K_OO,const Matrix& K_OI,const Matrix& K_IO,const Matrix& K_II,
//                       const Matrix& I_O,const Matrix& I_I)
//{
//    matE = V_OO-V_OI*V_II_fact.solve(V_IO);
//}
//
//void Domain::compute_tmp(const Matrix& V_OO,const Matrix& V_OI,const Matrix& V_IO,const Matrix& V_II,
//                       const Matrix& K_OO,const Matrix& K_OI,const Matrix& K_IO,const Matrix& K_II,
//                       const Matrix& I_O,const Matrix& I_I)
//{
//    mattmp = K_OO+0.5*I_O-V_OI*V_II_fact.solve(K_IO);
//}

void Domain::compute_helper_matrices()
{        
    Matrix V = SingleLayerOp();
    Matrix K = DoubleLayerOp();

    V_OO = V.block(0,0,N_outer_dofs,N_outer_dofs);
    V_OI = V.block(0,N_outer_dofs,N_outer_dofs,N_inner_dofs);
    V_IO = V.block(N_outer_dofs,0,N_inner_dofs,N_outer_dofs);
    V_II = V.block(N_outer_dofs,N_outer_dofs,N_inner_dofs,N_inner_dofs);
    K_OO = K.block(0,0,N_outer_dofs,N_outer_dofs);
    K_OI = K.block(0,N_outer_dofs,N_outer_dofs,N_inner_dofs);
    K_IO = K.block(N_outer_dofs,0,N_inner_dofs,N_outer_dofs);
    K_II = K.block(N_outer_dofs,N_outer_dofs,N_inner_dofs,N_inner_dofs);
    I_O = Matrix::Identity(N_outer_dofs,N_outer_dofs);
    I_I = Matrix::Identity(N_inner_dofs,N_inner_dofs);

    V_II_fact.compute(V_II);
    C_fact.compute(K_OO+0.5*I_O-V_OI*V_II_fact.solve(K_IO));
    matD = V_OI*V_II_fact.solve(K_II+0.5*I_I)-K_OI;
}

Vector Domain::compute_inner_deriv(const Vector& tracei, const Vector& derivo) const
{
    if(N_inner_bnds==0)
    {
        cerr<<"WARNING: no inner boundaries"<<endl;
        return derivo;
    }
   
    Vector Dtracei = V_OI*V_II_fact.solve((K_II+0.5*I_I)*tracei)-K_OI*tracei;
    Vector derivo_wm = add_outer_measure(derivo);
    Vector Ederivo = V_OO*derivo_wm-V_OI*V_II_fact.solve(V_IO*derivo_wm);

    Vector Atracei = V_II_fact.solve((K_II+0.5*I_I)*tracei+K_IO*C_fact.solve(Dtracei));
    Vector Bderivo = V_II_fact.solve(K_IO*C_fact.solve(Ederivo)-V_IO*derivo_wm);

    return remove_inner_measure(Atracei+Bderivo);        
}

Vector Domain::compute_outer_trace(const Vector& tracei, const Vector& derivo) const
{
    if(N_inner_bnds==0)
    {
        cerr<<"WARNING: no inner boundaries"<<endl;
        return tracei;
    }
  
//    Vector Dtracei = V_OI*V_II_fact.solve((K_II+0.5*I_I)*tracei)-K_OI*tracei;
    Vector derivo_wm = add_outer_measure(derivo);
    Vector Ederivo = V_OO*derivo_wm-V_OI*V_II_fact.solve(V_IO*derivo_wm);
//    Matrix D = V_OI*V_II_fact.solve(K_II+0.5*I_I)-K_OI;
//    Matrix E = V_OO-V_OI*V_II_fact.solve(V_IO);

    return C_fact.solve(matD*tracei+Ederivo);
    
}

Vector Domain::DirichletToNeumann(const Vector& trace) const
{
    return remove_measure(PoincareSteklow()*trace);
}

Vector Domain::NeumannToDirichlet(const Vector& deriv) const
{
    Vector trace = PoincareSteklow().ldlt().solve(add_measure(deriv));
    return trace.array()-trace.array().mean();
}

Real Domain::L2norm(const Vector& v) const
{
    Real norm=0.;
    for(unsigned bnd=0;bnd<N_bnds;bnd++)
        norm += pow(L2norm(v.segment(bnd_dofs[bnd],bnd_size[bnd]),bnd),2);

    return sqrt(norm);
}

Real Domain::L2norm_inner(const Vector& v) const
{
    Real norm=0.;
    unsigned bnd0 = N_bnds>1 ? 1:0;
    for(unsigned bnd=bnd0;bnd<N_bnds;bnd++)
        norm += pow(bnds[bnd].L2norm(v.segment(bnd_dofs[bnd]-bnd_dofs[bnd0],bnd_size[bnd])),2);
    return sqrt(norm);
}

Real Domain::L2norm(const Vector& v, unsigned i) const
{
    return bnds[i].L2norm(v);
}

Real Domain::Linfnorm(const Vector& v) const
{
    return v.lpNorm<Eigen::Infinity>();
}

Vector Domain::interpolate(const Vector& v, const Domain& dom) const
{
    Vector vi(N_outer_dofs+N_inner_dofs);
    
    for(unsigned i=0;i<N_bnds;i++)
        vi.segment(bnd_dofs[i],bnd_size[i]) = 
                bnds[i].interpolate(v.segment(dom.bnd_dofs[i],dom.bnd_size[i]),dom.bnds[i]);
    return vi;
}

Vector Domain::interpolate_inner(const Vector& v, const Domain& dom) const
{
    Vector vi(N_inner_dofs);
    unsigned bnd0 = N_bnds>1 ? 1:0;
    for(unsigned i=bnd0;i<N_bnds;i++)
        vi.segment(bnd_dofs[i]-bnd_dofs[bnd0],bnd_size[i]) = 
                bnds[i].interpolate(v.segment(dom.bnd_dofs[i]-dom.bnd_dofs[bnd0],dom.bnd_size[i]),dom.bnds[i]);
    return vi;
}

Vector Domain::remove_measure(const Vector& deriv) const
{
    Vector new_deriv(deriv.size());
    for(unsigned bnd=0;bnd<N_bnds;bnd++)
    {
        new_deriv.segment(bnd_dofs[bnd],bnd_size[bnd])=
            deriv.segment(bnd_dofs[bnd],bnd_size[bnd])
             .cwiseQuotient(bnds[bnd].get_measure());
    }
    return new_deriv;
}

Vector Domain::add_measure(const Vector& deriv) const
{
    Vector new_deriv(deriv.size());
    for(unsigned bnd=0;bnd<N_bnds;bnd++)
    {
        new_deriv.segment(bnd_dofs[bnd],bnd_size[bnd])=
            deriv.segment(bnd_dofs[bnd],bnd_size[bnd])
             .cwiseProduct(bnds[bnd].get_measure());
    }
    return new_deriv;
}

Vector Domain::remove_measure(const Vector& deriv, unsigned bnd) const
{   
    return deriv.cwiseQuotient(bnds[bnd].get_measure());
}

Vector Domain::add_measure(const Vector& deriv, unsigned bnd) const
{
    return deriv.cwiseProduct(bnds[bnd].get_measure());
}

Vector Domain::remove_inner_measure(const Vector& deriv) const
{
    Vector new_deriv(deriv.size());
    unsigned bnd0 = N_bnds>1? 1:0;
    for(unsigned bnd=bnd0;bnd<N_bnds;bnd++)
    {
        new_deriv.segment(bnd_dofs[bnd]-bnd_dofs[bnd0],bnd_size[bnd])=
            deriv.segment(bnd_dofs[bnd]-bnd_dofs[bnd0],bnd_size[bnd])
             .cwiseQuotient(bnds[bnd].get_measure());
    }
    return new_deriv;
}

Vector Domain::add_inner_measure(const Vector& deriv) const
{
    Vector new_deriv(deriv.size());
    unsigned bnd0 = N_bnds>1? 1:0;
    for(unsigned bnd=bnd0;bnd<N_bnds;bnd++)
    {
        new_deriv.segment(bnd_dofs[bnd]-bnd_dofs[bnd0],bnd_size[bnd])=
            deriv.segment(bnd_dofs[bnd]-bnd_dofs[bnd0],bnd_size[bnd])
             .cwiseProduct(bnds[bnd].get_measure());
    }
    return new_deriv;
}

Vector Domain::remove_outer_measure(const Vector& deriv) const
{   
    return deriv.cwiseQuotient(bnds[0].get_measure());
}

Vector Domain::add_outer_measure(const Vector& deriv) const
{
    return deriv.cwiseProduct(bnds[0].get_measure());
}

Vector Domain::get_outer(const Vector& v) const
{
    return v.head(N_outer_dofs);
}

Vector Domain::get_inner(const Vector& v) const
{
    return v.tail(N_inner_dofs);
}

void Domain::mesh()
{
    if(meshed)
        return;
    
    gmsh::initialize();
    
    string model_name = "dom";
    
    gmsh::model::add(model_name);
    
    vector<vector<int>> pts_tags(N_bnds);
    vector<int> bnd_tags(N_bnds);
    vector<int> curve_loops_tags(N_bnds);

    Real point_mesh_size = 1.*h_min; // to add a mesh constraint, here below use addPoint(G(i,0),G(i,1),0.,point_mesh_size)
    for(auto bnd = 0u; bnd<N_bnds;bnd++)
    {
        const Matrix& G = bnds[bnd].get_G();
        unsigned n_pts = bnds[bnd].get_n_pts();
        for(auto i=0u;i<n_pts;i++)
        {
            point_mesh_size = (G.row((i-1+n_pts)%n_pts)-G.row((i+1)%n_pts)).norm();
            pts_tags[bnd].push_back(gmsh::model::geo::addPoint(G(i,0),G(i,1),0.,point_mesh_size));
        }
        pts_tags[bnd].push_back(pts_tags[bnd].front());

        bnd_tags[bnd] = gmsh::model::geo::addPolyline(pts_tags[bnd]);
        curve_loops_tags[bnd] = gmsh::model::geo::addCurveLoop({bnd_tags[bnd]});
    }
        
    unsigned domain_surface= gmsh::model::geo::addPlaneSurface(curve_loops_tags);        
    
    gmsh::model::geo::synchronize();
    
    gmsh::model::mesh::generate(2);
    
//    gmsh::write(model_name+".msh");
    
    
    // now extract points and elements from the generated mesh
    gmsh::vectorpair ents;
    gmsh::model::getEntities(ents);
    
    vector<vector<size_t>> elems_tags;
    size_t max_tag=0;
    for(auto ent:ents)
    {
        int dim=ent.first;
        int tag=ent.second;
        
        if(dim<2)
            continue;
        
        std::vector<int> elemTypes;
        std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);
        
        for(auto elemT:elemNodeTags) 
        for(auto i=0u;i<elemT.size()/3;i++)
        {
            max_tag = max(max_tag,elemT[3*i]);
            max_tag = max(max_tag,elemT[3*i+1]);
            max_tag = max(max_tag,elemT[3*i+2]);
            elems_tags.push_back({elemT[3*i],elemT[3*i+1],elemT[3*i+2]});
        }
    }
    
    vector<bool> used_nodes(max_tag,false);
    for(auto elem:elems_tags)
    for(auto node:elem)
        used_nodes[node-1]=true;
    
    unsigned n_nodes=0;
    for(auto node:used_nodes)
        if(node)
            n_nodes++;
    
    points.resize(2*n_nodes);
    vector<int> tag_to_ind(max_tag,-1);
    unsigned ind=0;
    for(auto ent:ents)
    {
        int dim=ent.first;
        int tag=ent.second;
        
        std::vector<std::size_t> nodeTags;
        std::vector<double> nodeCoords, nodeParams;
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);
        
        for(auto i=0u;i<nodeTags.size();i++)
        {
            if(used_nodes[nodeTags[i]-1])
            {
                tag_to_ind[nodeTags[i]-1]=ind;
                points[2*ind]=nodeCoords[3*i];
                points[2*ind+1]=nodeCoords[3*i+1];
                ind++;
            }
        }
    }
    
    elems.clear();
    for(auto elem:elems_tags)
    for(auto node:elem)
        elems.push_back(tag_to_ind[node-1]);
    
    gmsh::finalize();
    
    //detect points close to boundary. Evaluating the s/l op at these points is dangerous
    is_point_on_bnd = vector<int>(points.size()/2,-1);
    map_to_bnd_points.clear();
    Real x,y;
    unsigned counter_bnd_pts=0;
    bool is_on_bnd;
    Point p,q,r;
    Real t,s;
    for(auto i=0u;i<is_point_on_bnd.size();i++)
    {
        is_on_bnd=false;
        
        x = points[2*i];
        y = points[2*i+1];
        r = {x,y};
        
        for(auto bnd=0u;bnd<N_bnds;bnd++)
        {
            unsigned n_pts = bnds[bnd].get_n_pts();
            const Matrix& G = bnds[bnd].get_G();
            p = G.bottomRows<1>();

            for(auto j=0u;j<G.rows();j++)
            {
                q = G.row(j);
                if(is_point_on_segment(p,q,r,t,s))
                {
                    is_point_on_bnd[i] = counter_bnd_pts;
                    unsigned indq=j;
                    unsigned indp=(j+n_pts-1)%n_pts;

                    map_to_bnd_points.push_back({bnd,indp,indq,t,s});
                    counter_bnd_pts++;

                    is_on_bnd=true;
                    break;
                }
                p=q;
            }
            if(is_on_bnd)
                break;
        }        
    }
    
    boundary_pts.clear();
    non_boundary_pts.clear();
    boundary_pts.reserve(is_point_on_bnd.size());
    non_boundary_pts.reserve(is_point_on_bnd.size());
    for(unsigned i=0;i<is_point_on_bnd.size();i++)
    {
        if(is_point_on_bnd[i]>=0)
            boundary_pts.push_back(i);
        else
            non_boundary_pts.push_back(i);
    }
    boundary_pts.shrink_to_fit();
    non_boundary_pts.shrink_to_fit();
    
    meshed=true;
 }

bool Domain::is_point_on_segment(const Point& p, const Point& q, 
                                 const Point& r, Real& t, Real& s)
{
    Real tol = 1.;//1e-2;
    Point qmp = q-p;
    Point rmp = r-p;
    Point n{qmp(1),-qmp(0)};
    Real norm_n = n.norm();
    t = (rmp(0)*qmp(0)+rmp(1)*qmp(1))/norm_n/norm_n;
    
    if(0<=t && t<=1)
    {
        s = (n(0)*rmp(0)+n(1)*rmp(1))/norm_n/norm_n;
        Real dist = abs(s*norm_n);
        
        if(dist<tol*norm_n)
            return true;
        else
            return false;
    }
    
    return false;
}

void Domain::write_sol_and_mesh(const string sol_name, const vector<Real>& u) const
{
    if(points.empty() || elems.empty())
    {
        cerr<<"Error: call Domain::mesh() first."<<endl;
        return;
    }
    
    // write a VTK solution
    const int dim = 2;
    const int cell_size = 3;
   
    string extension = ".vtu";
    string file_name = sol_name + extension;

    leanvtk::VTUWriter writer;
    writer.add_scalar_field("U", u);
    writer.write_surface_mesh(file_name, dim, cell_size, points, elems);
}

Real Domain::eval_on_point(const Point p, const Vector& trace, const Vector& deriv) const
{
    // WARNING:  we suppose that p doesn't lie on a boundary and therefore
    // dont mind about singularities.
    
    Eigen::ArrayXd diffx, diffy;
    Eigen::ArrayXd V,K,normsq;
    
    Real u_num = 0.;
    
    for(auto b=0u;b<N_bnds;b++)
    {
        const Matrix& G = bnds[b].get_G();
        const Matrix& N = bnds[b].get_N();
        const Vector& measure = bnds[b].get_measure();
        const unsigned n = bnds[b].get_n_pts();
        
//        for(unsigned i=0;i<n;i++)
//        {
//            Real V,K;
//            Point q,Ni;
//            Real derivi,tracei,measurei;
//            
//            q = 0.5*(G.row(i)+G.row((i+1)%n));
//            Ni = 0.5*(N.row(i)+N.row((i+1)%n));
//            measurei = 0.5*(measure(i)+measure((i+1)%n));
//            derivi = 0.5*(deriv[bnd_dofs[b]+i]+deriv[bnd_dofs[b]+(i+1)%n]);
//            tracei = 0.5*(trace[bnd_dofs[b]+i]+trace[bnd_dofs[b]+(i+1)%n]);
//            
////            q = G.row(i);
////            Ni = N.row(i);
////            measurei = measure(i);
////            derivi = deriv[bnd_dofs[b]+i];
////            tracei = trace[bnd_dofs[b]+i];
//            
//            V = log((p-q).squaredNorm());
//            K = ((p(0)-q(0))*Ni(0)+(p(1)-q(1))*Ni(1))/(p-q).squaredNorm();
//            
//            V = (-1./4./M_PI)*V*derivi*measurei;
//            K = (-1./2./M_PI)*K*tracei*measurei;
//            
//            u_num += (V+K);
//        }
//        
        //vectorized version
        diffx = p(0)-G.col(0).array();
        diffy = p(1)-G.col(1).array();
        normsq = diffx.square()+diffy.square();
        V = normsq.log();
        K = (diffx*N.col(0).array()+diffy*N.col(1).array())/normsq;
        V = (-1./4./M_PI)*V*deriv.segment(bnd_dofs[b],bnd_size[b]).array()*measure.array();
        K = (-1./2./M_PI)*K*trace.segment(bnd_dofs[b],bnd_size[b]).array()*measure.array();
        u_num += V.sum()+K.sum();
    }
    
    return u_num;
}

vector<Real> Domain::eval_on_points(const vector<Real>& pts, const Vector& trace, const Vector& deriv) const
{
    unsigned n_points = pts.size()/2;
    vector<Real> u(n_points,0.);
    
    for(auto p=0u;p<n_points;p++)
        u[p] = eval_on_point({pts[2*p],pts[2*p+1]},trace,deriv);
     
    return u;
}

vector<Real> Domain::get_non_boundary_points() const
{
    if(points.empty() || elems.empty())
    {
        cerr<<"Error: call Domain::mesh() first."<<endl;
        return vector<Real>();
    }
    
    unsigned n_points = points.size()/2;
    vector<Real> pts;
    
    for(auto p=0u;p<n_points;p++)
    {        
        if(is_point_on_bnd[p]<0) //then x,y isnt on boundary
        {
            pts.push_back(points[2*p]);
            pts.push_back(points[2*p+1]);
        }        
    }
    
    return pts;
}

vector<Real> Domain::eval_on_mesh(const Vector& trace, const Vector& deriv) const
{
    if(points.empty() || elems.empty())
    {
        cerr<<"Error: call Domain::mesh() first."<<endl;
        return vector<Real>();
    }
    
    unsigned n_points = points.size()/2;
    vector<Real> u(n_points,0.);
    
//    #pragma omp parallel for
//    for(auto p=0u;p<n_points;p++)
//    {        
//        if(is_point_on_bnd[p]>=0) //then x,y is on boundary
//        {
//            const tuple<unsigned,unsigned,unsigned,Real,Real>& to_bnd_pts = map_to_bnd_points[is_point_on_bnd[p]];
//            const unsigned bnd = std::get<0>(to_bnd_pts);
//            const unsigned indp = std::get<1>(to_bnd_pts);
//            const unsigned indq = std::get<2>(to_bnd_pts);
//            const Real t = std::get<3>(to_bnd_pts);
//            const Real s = std::get<4>(to_bnd_pts);
//                        
//            Real trace_p = trace[bnd_dofs[bnd]+indp];
//            Real trace_q = trace[bnd_dofs[bnd]+indq];
//            u[p] = trace_p+t*(trace_q-trace_p);
//            
//            Real deriv_p = deriv[bnd_dofs[bnd]+indp];
//            Real deriv_q = deriv[bnd_dofs[bnd]+indq];
//            Real norm = (bnds[bnd].get_G(indp)-bnds[bnd].get_G(indq)).norm();
//            if(bnd==0)
//                u[p] += (deriv_p+t*(deriv_q-deriv_p))*s*norm;
//            else
//                u[p] -= (deriv_p+t*(deriv_q-deriv_p))*s*norm;
//        }
//        else
//            u[p] = eval_on_point({points[2*p],points[2*p+1]},trace,deriv);
//    }
    
    unsigned i,p;
    
    #pragma omp parallel for schedule(dynamic) default(shared) private(i,p)
    for(auto i=0u;i<boundary_pts.size();i++)
    {       
        p = boundary_pts[i];
        
        const tuple<unsigned,unsigned,unsigned,Real,Real>& to_bnd_pts = map_to_bnd_points[is_point_on_bnd[p]];
        const unsigned bnd = std::get<0>(to_bnd_pts);
        const unsigned indp = std::get<1>(to_bnd_pts);
        const unsigned indq = std::get<2>(to_bnd_pts);
        const Real t = std::get<3>(to_bnd_pts);
        const Real s = std::get<4>(to_bnd_pts);

        Real trace_p = trace[bnd_dofs[bnd]+indp];
        Real trace_q = trace[bnd_dofs[bnd]+indq];
        u[p] = trace_p+t*(trace_q-trace_p);

        Real deriv_p = deriv[bnd_dofs[bnd]+indp];
        Real deriv_q = deriv[bnd_dofs[bnd]+indq];
        Real norm = (bnds[bnd].get_G(indp)-bnds[bnd].get_G(indq)).norm();
        if(bnd==0)
            u[p] += (deriv_p+t*(deriv_q-deriv_p))*s*norm;
        else
            u[p] -= (deriv_p+t*(deriv_q-deriv_p))*s*norm;
        
    }
    
    #pragma omp parallel for schedule(dynamic) default(shared) private(i,p)
    for(i=0u;i<non_boundary_pts.size();i++)
    {
        p = non_boundary_pts[i];
        u[p] = eval_on_point({points[2*p],points[2*p+1]},trace,deriv);
    }
    
    return u;
}

ostream& Domain::write(ostream& os, string name) const
{
    
    os<<name<<" = cell("<<N_bnds<<",1);"<<endl;
    for(auto i=0u;i<N_bnds;i++)
        bnds[i].to_matlab(os,name+"{"+to_string(i+1)+"}");
    
    return os;
}

ostream& Domain::write_boundary_sol_and_coords(ostream& os, const Vector& v, string name) const
{
    Matrix xyu(N_outer_dofs+N_inner_dofs,3);
    for(unsigned i=0;i<N_bnds;i++)
    {
        xyu.block(bnd_dofs[i],0,bnd_size[i],2) = bnds[i].get_G();
        xyu.block(bnd_dofs[i],2,bnd_size[i],1) = v.segment(bnd_dofs[i],bnd_size[i]);
    }
    write_matrix(os,xyu,name);
    return os;
}

Point Domain::point(unsigned ind) const
{
    return {points[2*ind],points[2*ind+1]};
}

const vector<Boundary>& Domain::get_bnds() const
{
    return bnds;
}

Real Domain::get_h_min() const
{
    return h_min;
}

Real Domain::get_h_max() const
{
    return h_max;
}

unsigned Domain::get_tot_n_dofs() const
{
    return N_outer_dofs+N_inner_dofs;
}

void Domain::save(string name) const
{
    ofstream ofile(name+".bin", std::ios::out | std::ios::binary | std::ios::trunc);
    save(ofile);
    ofile.close();
}

void Domain::load(string name)
{
    std::ifstream ifile(name+".bin", std::ios::in | std::ios::binary);
    load(ifile);
    ifile.close();
}

void Domain::save(ofstream& out) const
{
    write_scalar_binary(out, N_bnds);
    write_scalar_binary(out, N_inner_bnds);
    write_scalar_binary(out, N_outer_dofs);
    write_scalar_binary(out, N_inner_dofs);
    write_scalar_binary(out, h_min);
    write_scalar_binary(out, h_max);
    
    for(unsigned i=0;i<N_bnds;i++)
        bnds[i].save(out);
    
    write_binary(out,bnd_dofs);
    write_binary(out,bnd_size);
    
    write_binary(out,points);
    write_binary(out,elems);
    write_binary(out,is_point_on_bnd);
    write_binary(out,map_to_bnd_points);
}

void Domain::load(ifstream& in)
{
    read_scalar_binary(in, N_bnds);
    read_scalar_binary(in, N_inner_bnds);
    read_scalar_binary(in, N_outer_dofs);
    read_scalar_binary(in, N_inner_dofs);
    read_scalar_binary(in, h_min);
    read_scalar_binary(in, h_max);
    
    bnds.clear();
    for(unsigned i=0;i<N_bnds;i++)
        bnds.push_back(Boundary(in));
    
    read_binary(in,bnd_dofs);
    read_binary(in,bnd_size);
    
    read_binary(in,points);
    read_binary(in,elems);
    read_binary(in,is_point_on_bnd);
    read_binary(in,map_to_bnd_points);
}

bool Domain::operator==(const Domain& dom) const
{
    if(N_bnds!=dom.N_bnds)
        return false;
    if(N_inner_bnds!=dom.N_inner_bnds)
        return false;
    if(N_outer_dofs!=dom.N_outer_dofs)
        return false;
    if(N_inner_dofs!=dom.N_inner_dofs)
        return false;
    if(h_min!=dom.h_min)
        return false;
    if(h_max!=dom.h_max)
        return false;
    
    for(unsigned i=0;i<N_bnds;i++)
        if(bnds[i]!=dom.bnds[i])
            return false;
    
    if(bnd_dofs!=dom.bnd_dofs)
        return false;
    if(bnd_size!=dom.bnd_size)
        return false;
    
    if(points!=dom.points)
        return false;
    if(elems!=dom.elems)
        return false;
    if(is_point_on_bnd!=dom.is_point_on_bnd)
        return false;
    if(map_to_bnd_points!=dom.map_to_bnd_points)
        return false;
    
    return true;
}

bool Domain::operator%(const Domain& dom) const
{
    if(N_bnds!=dom.N_bnds)
        return false;
    if(N_inner_bnds!=dom.N_inner_bnds)
        return false;
    if(N_outer_dofs!=dom.N_outer_dofs)
        return false;
    if(N_inner_dofs!=dom.N_inner_dofs)
        return false;
    
    Real tol = 1e-8;
    if(abs(h_min-dom.h_min)>tol*min(h_min,dom.h_min))
        return false;
    if(abs(h_max-dom.h_max)>tol*min(h_max,dom.h_max))
        return false;
    
    if(bnd_dofs!=dom.bnd_dofs)
        return false;
    if(bnd_size!=dom.bnd_size)
        return false;
    
    for(unsigned i=0;i<N_bnds;i++)
        if(!(bnds[i]%dom.bnds[i]))
            return false;
    
    return true;
}

ostream& operator<<(ostream& os, const uBidomain::Domain& dom)
{
    for(auto bnd=0u;bnd<dom.get_bnds().size();bnd++)
        os<<dom.get_bnds()[bnd]<<endl;
    return os;
}