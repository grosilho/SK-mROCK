#include <limits>

#include "Boundary.h"
#include "Utils.h"

using namespace uBidomain;

Boundary::Boundary()
{
    
}

Boundary::Boundary(vector<Side> sides)
{
    smooth_boundary = false;
    
    n_pts=0;
    for(auto side:sides)
        n_pts += side.get_n();
    
    t = Eigen::ArrayXd::LinSpaced(n_pts, 0, n_pts - 1) / n_pts;
    dt = (t(1)-t(0))*Vector::Ones(n_pts);
    t = t+dt/2.;
    
    G = Matrix::Zero(n_pts,2);
    dG = Matrix::Zero(n_pts,2);
    ddG = Matrix::Zero(n_pts,2);
    N = Matrix::Zero(n_pts,2);
    norm_dG = Vector::Zero(n_pts);
    norm_ddG = Vector::Zero(n_pts);
    h = std::numeric_limits<Real>::max();
    length = 0.;
    
    unsigned n_sides = sides.size();
    unsigned dofs_count=0;
    
    Real tol=1e-10;
    Point q = sides.back().end();
    int prev_side = -1;
    int current_side;
    sides_dofs.resize(n_sides);
    Real dt_sc = dt(0);
    
    for(auto i=0u;i<n_sides;i++)
    {
        //search for the side whose starting point matches the previous side ending point
        for(int j=prev_side+1;j<prev_side+1+n_sides;j++)
        {
            current_side = j%n_sides;
            if((sides[current_side].begin()-q).norm()<tol)
                break;
        }
        const Side& side = sides[current_side];
        prev_side=current_side;
        q = side.end();
        
        unsigned side_size = side.get_n();
        
        sides_dofs[i][0]=dofs_count;
        sides_dofs[i][1]=side_size;
        sides_dofs[i][2]=current_side; //stores the original index
        
        Real dx = side.get_dx();
        G.block(dofs_count,0,side_size,2) = side.get_G();
        dG.block(dofs_count,0,side_size,2) = (dx/dt_sc)*side.get_dG();
        ddG.block(dofs_count,0,side_size,2) = (dx/dt_sc)*(dx/dt_sc)*side.get_ddG();
        N.block(dofs_count,0,side_size,2) = side.get_N();
        norm_dG.segment(dofs_count,side_size) = (dx/dt_sc)*side.get_norm_dG();
        norm_ddG.segment(dofs_count,side_size) = (dx/dt_sc)*(dx/dt_sc)*side.get_norm_ddG();
        h = min(h,side.get_h());
        length += side.get_length();
                        
        dofs_count += side_size;
    }
    
    if(smooth_boundary)
        fourier_interpolation();

    measure = norm_dG.array()*dt.array();

    //detect corners by comparing normals.
    Real angle_tol = 0.3;
    for(auto side=0u;side<n_sides;side++)
    {
        //check beginning of this side
        unsigned ind1 = sides_dofs[side][0];
        unsigned ind2 = (ind1-1+n_pts)%n_pts;
        Real dotp = N(ind1,0)*N(ind2,0)+N(ind1,1)*N(ind2,1);
        Real alpha = acos(dotp);

        if(alpha>angle_tol)
            corners.push_back({ind1,ind2});

        //check end of this side
        ind1 = sides_dofs[side][0]+sides_dofs[side][1]-1;
        ind2 = (ind1+1+n_pts)%n_pts;
        dotp = N(ind1,0)*N(ind2,0)+N(ind1,1)*N(ind2,1);
        alpha = acos(dotp);

        if(alpha>angle_tol)
            corners.push_back({ind1,ind2});
    }
    
}

Boundary::Boundary(string filename)
{
    load(filename);
}

Boundary::Boundary(ifstream& in)
{
    load(in);
}

Boundary Boundary::get_circle(Point c, Real r, Real dG, Real sign_N, Real kappa)
{
    return Boundary({Arc(c,r,0.,2*M_PI,dG,sign_N,kappa)});
}

Boundary Boundary::get_ellipse(Point c, Real r1, Real r2, Real dG, Real sign_N, Real kappa)
{
    return Boundary({EllipseArc(c,r1,r2,0.,2*M_PI,dG,sign_N,kappa)});
}

Boundary Boundary::get_square(Point c, Real s, Real dG, Real sign_N, Real kappa)

{
    return Boundary({c+(s/2.)*Segment({1.,-1.},{1.,1.},2.*dG/s,sign_N,kappa),
                     c+(s/2.)*Segment({1.,1.},{-1.,1.},2.*dG/s,sign_N,kappa),
                     c+(s/2.)*Segment({-1.,1.},{-1.,-1.},2.*dG/s,sign_N,kappa),
                     c+(s/2.)*Segment({-1.,-1.},{1.,-1.},2.*dG/s,sign_N,kappa)});
}

Boundary Boundary::get_rectangle(Point p1, Point p2, Real dG, Real sign_N, Real kappa)

{
    Point p3(p2(0),p1(1));
    Point p4(p1(0),p2(1));
    return Boundary({Segment(p1,p3,dG,sign_N,kappa),
                     Segment(p3,p2,dG,sign_N,kappa),
                     Segment(p2,p4,dG,sign_N,kappa),
                     Segment(p4,p1,dG,sign_N,kappa)});
}

void Boundary::find_closest_point(const Point& p, Real& t0, Real& err) const
{    
    Real tolN = 1e-14;
    Point q,dq,ddq;
    Real c, dc;
    
    q(0) = FGx.eval_one_pt(t0);
    q(1) = FGy.eval_one_pt(t0);   
                      
    err = (p-q).norm();
    Real dtq = 1.;
    
    unsigned iter=0;
    
//    cout<<"t = "<<t0<<endl;
//    cout<<"p = "<<p<<endl;
//    cout<<"q = "<<q<<endl;
//    cout<<"err = "<<err<<endl;        

    while(err>tolN*(1.+p.norm()) && abs(dtq)>tolN)
    {
        iter++;
        
        dq(0) = dFGx.eval_one_pt(t0);
        dq(1) = dFGy.eval_one_pt(t0);
        ddq(0) = ddFGx.eval_one_pt(t0);
        ddq(1) = ddFGy.eval_one_pt(t0);
        
        c = dq.dot(p-q);
        dc = -dq.squaredNorm()+ddq.dot(p-q);
        dtq = -c/dc;

        t0 += dtq;
        q(0) = FGx.eval_one_pt(t0);
        q(1) = FGy.eval_one_pt(t0);
        err = (p-q).norm();
        
//        cout<<"dtp = "<<dtq<<endl;
//        cout<<"q = "<<q<<endl;
//        cout<<"dq = "<<dq<<endl;
//        cout<<"ddq = "<<ddq<<endl;
//        cout<<"Err = "<<err<<endl;
//        int kk;
//        if(iter>10)
//            cin>>kk;
    }
  
}

void Boundary::fourier_interpolation()
{
    //as it is, the boundary is usually non smooth if it is constructed gluying 
    //different Side. This is due to corners and discontinuity of norm_dG
    //the purpose of this method is to smooth everything out
    //to do so we take the collocation nodes G, interpolate them with Fourier
    //and redefine normals and derivatives
    
    Real tolF = 1e-8; //maximal Fourier interpolation error
    Real N_min = 1e2;//minimal number of Fourier modes
    
    Matrix tmpG(G.rows(),G.cols());
    unsigned dn = G.rows()/2;
    dn = min(1100u,dn);
    tmpG.topRows(dn) = G.bottomRows(dn);
    tmpG.bottomRows(G.rows()-dn) = G.topRows(G.rows()-dn);
    Matrix fc = computeFcoeff(tmpG, tolF, N_min);
    
//    unsigned dn=0;
//    Matrix fc = computeFcoeff(G, tolF, N_min);

    FGx.set_FCoeff(fc.col(0));
    FGy.set_FCoeff(fc.col(1));
        
    dFGx.set_FCoeff(fc.col(0));
    dFGy.set_FCoeff(fc.col(1));
    dFGx.diff();
    dFGy.diff();
        
    ddFGx.set_FCoeff(fc.col(0));
    ddFGy.set_FCoeff(fc.col(1));
    ddFGx.diff().diff();
    ddFGy.diff().diff();    

    //for every point G find the closes point on the interpolated curve    
    Real err;
    Real last_dt=t(1)-t(0);
    Real t0 = dn*last_dt;
//    t = t.array()+t0;
    Real max_err = 0.;   
    for(unsigned i=0;i<n_pts;i++)
    {
        find_closest_point(G.row(i),t0,err); 
        t(i) = t0;                                
        if(i>0)
            last_dt = t(i)-t(i-1);
        t0 += last_dt;
        
//        find_closest_point(G.row(i),t(i),err);
        max_err = max(max_err,err);
    }

//    cout<<"Max interpolation error = "<<max_err<<endl;
//    cout<<"norm t "<<t.norm()<<endl;
//    bool linear_t=false; 
//    if(linear_t)
//        t = Vector::LinSpaced(n_pts,0,n_pts-1.)/n_pts;     
//    else
//    {
        for(unsigned i=1;i<t.size();i++)
            if(t(i)<=t(i-1))
                cerr<<"ERROR: in fourier interpolation, new t is not increasing"<<endl;
        if(abs(t(0)-t(t.size()-1))>=1.)
            cerr<<"ERROR: in fourier interpolation, t has range greater than 1"<<endl;
//    }   
        
    dt.resize(n_pts);
    for(unsigned i=1;i<dt.size()-1;i++)
        dt(i) = (t(i+1)-t(i-1))/2.;
    dt(0) = ( t(1)-t(n_pts-1)+1. )/2.;
    dt(n_pts-1) = (t(0)+1.-t(n_pts-2))/2.;
        
    G.col(0) = FGx.eval(t);
    G.col(1) = FGy.eval(t);
    
    // get the sign of normal
    Point N0 = N.row(0);
    Point N1;
    N1(0) = dG(0,1)/norm_dG(0);
    N1(1) = -dG(0,0)/norm_dG(0);
    Real cosN0N1 = N0.dot(N1);
    Real signN;
    if(abs(cosN0N1-1.)<0.1)
        signN = 1.;
    else if(abs(cosN0N1+1.)<0.1)
        signN = -1.;
    else
        cout<<"ERROR in Boundary fourier_interpolation"<<endl;
        
    dG.col(0) = dFGx.eval(t);
    dG.col(1) = dFGy.eval(t);
    
    norm_dG = (dG.col(0).array().pow(2)+dG.col(1).array().pow(2)).sqrt();

    N.col(0) = signN*dG.col(1).array()/norm_dG.array();
    N.col(1) = -signN*dG.col(0).array()/norm_dG.array();
       
    ddG.col(0) = ddFGx.eval(t);
    ddG.col(1) = ddFGy.eval(t);
    
    norm_ddG = (ddG.col(0).array().pow(2)+ddG.col(1).array().pow(2)).sqrt();
    
//    cout<<"Old length = "<<length<<endl;
//    cout<<"New length = "<<(norm_dG.array()*new_dt.array()).sum()<<endl;
    
}

Matrix Boundary::computeFcoeff(Matrix& Gamma, Real tol, int Nmin)
{
    Vector fcx = FourierInterpolation(Gamma.block(0,0,Gamma.rows(),1), tol, Nmin);
    Vector fcy = FourierInterpolation(Gamma.block(0,1,Gamma.rows(),1), tol, Nmin);
    int dim = max(fcx.size(),fcy.size());
    Matrix fcG = Matrix::Zero(dim,2);
    fcG.block((dim-fcx.size())/2,0,fcx.size(),1) = fcx;
    fcG.block((dim-fcy.size())/2,1,fcy.size(),1) = fcy;
    return fcG;
}

Vector Boundary::FourierInterpolation(Vector curve, Real tol, int Nmin)
{
    FSeries fs;
    fs.interp(curve);
    fs.truncate(tol, Nmin);
    return fs.get_FCoeff();
}

const unsigned Boundary::get_n_pts() const
{
    return n_pts;
}

const Vector& Boundary::get_t() const
{
    return t;
}

const Vector& Boundary::get_dt() const
{
    return dt;
}

Real Boundary::get_h() const
{
    return h;
}

Real Boundary::get_length() const
{
    return length;
}

const Matrix& Boundary::get_G() const
{
    return G;
}

const Point Boundary::get_G(unsigned i) const
{
    return G.block(i,0,1,2);
}

const Matrix& Boundary::get_dG() const
{
    return dG;
}

const Matrix& Boundary::get_ddG() const
{
    return ddG;
}

const Matrix& Boundary::get_N() const
{
    return N;
}

const Point Boundary::get_N(unsigned i) const
{
    return N.block(i,0,1,2);
}

const Vector& Boundary::get_norm_dG() const
{
    return norm_dG;
}

const Vector& Boundary::get_norm_ddG() const
{
    return norm_ddG;
}

const Vector& Boundary::get_measure() const
{
    return measure;
}

const Real Boundary::get_measure(unsigned i) const
{
    return measure(i);
}

ostream& Boundary::to_matlab(ostream& os, string name) const
{
    if(!name.empty())
        name += ".";
    
    os<<name<<"n = "<<n_pts<<";"<<endl;
    
    write_matrix(os,t,name+"x");
    write_matrix(os,G,name+"G");
    write_matrix(os,dG,name+"dG");
    write_matrix(os,ddG,name+"ddG");
    write_matrix(os,N,name+"N");
    write_matrix(os,norm_dG,name+"norm_dG");
    write_matrix(os,norm_ddG,name+"norm_ddG");   
    
    return os;
}

vector<array<unsigned,2>> Boundary::get_corners() const
{
    return corners;
}

vector<array<unsigned,3>> Boundary::get_sides_dofs() const
{
    return sides_dofs;
}

Vector Boundary::interpolate(const Vector& v, const Boundary& bnd) const
{
    if(smooth_boundary)
        return interpolate_smooth(v,bnd);
    else
        return interpolate_nonsmooth(v,bnd);
}

void Boundary::find_closest_discretized_point(const Real t0, Real& tj, unsigned& j, const Vector& dtc) const
{
    //find  j=argmin |t(j)-t0| with t(j)<=t0
    //the argument j is a first guess
    
    Real tol=1e-13;
    while(tj>t0)
    {
        j = ((j+n_pts)-1)%n_pts;
        tj -= dtc(j);
    }
    while(tj+dtc(j)<t0+tol/2.)
    {        
        tj += dtc(j);
        j = (j+1)%n_pts;
    }
    
    //here we just remove the accumulated rounding errors
    if(abs(tj-t(j))<0.5)   
        tj=t(j);  
    else if(tj<t(j))
        tj = t(j)-1.;
    else if(tj>t(j))
        tj = t(j)+1.;
    
    if(abs(t(j)-tj)>tol && abs(abs(t(j)-tj)-1.)>tol)
    {
        cout<<"ERROR n Boundary::find_closest_discretized_point"<<endl;
//        cout<<"tj = "<<tj<<", t(j) = "<<t(j)<<", err = "<<abs(t(j)-tj)<<endl;
    }
}

Vector Boundary::interpolate_smooth(const Vector& v, const Boundary& bnd) const
{   
    // inteprolate v, living on bnd, onto "this"
    //to do so we use Newton divided difference formula
    const unsigned order = 4; //intepolation order
    
    Vector v_interp = Vector::Zero(n_pts);
    Real tol = 1e-13;    
    
    vector<unsigned> ind_j(order);
    vector<Real> t_j(order);
    vector<Real> dt_j(order);
    vector<vector<Real>> Ncoeff(order,vector<Real>(order));
    int sign;
    
    const Vector& tc = bnd.t;
    const Vector& tf = this->t;
    
    unsigned n_pts_c = bnd.n_pts;
    unsigned n_pts_f = this->n_pts;
    
    Vector dtc(n_pts_c);
    dtc.segment(0,n_pts_c-1) = tc.segment(1,n_pts_c-1)-tc.segment(0,n_pts_c-1);
    dtc(n_pts_c-1) = 1-tc(n_pts_c-1)+tc(0);
    
    Vector dtf(n_pts_f);
    dtf.segment(0,n_pts_f-1) = tf.segment(1,n_pts_f-1)-tf.segment(0,n_pts_f-1);
    dtf(n_pts_f-1) = 1-tf(n_pts_f-1)+tf(0);       
    
    Real s0,s1,t1;        
    s1 = tc(0)+ dtf(0)/2.-dtc(0)/2.;
    Real last_ds = 0.;
    Real err;
    unsigned j1=0;
    t1 = tc(j1);
    
    for(unsigned i=0;i<n_pts_f;i++)
    {
        s0 = s1;
        s1 += last_ds;  //a first guess
        bnd.find_closest_point(G.row(i),s1,err);
        bnd.find_closest_discretized_point(s1,t1,j1,dtc);
        
        if(abs(s1-t1)<tol)
            v_interp(i)=v(j1);
        else
        {                                                
            unsigned j2=(j1+1)%n_pts_c;
            Real t2 = t1 + dtc(j1);     
            
            if(abs(s1-t1)<abs(s1-t2))
            {
                ind_j[0] = j1;
                ind_j[1] = j2;
                t_j[0] = t1;
                t_j[1] = t2;
                sign = -1.;
            }
            else
            {
                ind_j[0] = j2;
                ind_j[1] = j1;
                t_j[0] = t2;
                t_j[1] = t1;
                sign = 1.;
            }
            for(unsigned o=2;o<order;o++)
            {
                ind_j[o] = (ind_j[o-2]+sign+n_pts_c)%n_pts_c;
                t_j[o] = t_j[o-2]+sign*dtc[ind_j[o-(sign+1)]];
                sign=-sign;
            }
            
            for(unsigned l=0;l<order;l++)
                Ncoeff[l][0] = v(ind_j[l]);
            
            for(unsigned c=1;c<order;c++)
            for(unsigned l=c;l<order;l++)
                Ncoeff[l][c] = (Ncoeff[l][c-1]-Ncoeff[l-1][c-1])
                                /(t_j[l]-t_j[l-c]);
                    
            v_interp(i) = Ncoeff[0][0];
            Real wk = 1.;
            for(unsigned k=1;k<order;k++)
            {
                wk *= (s1-t_j[k-1]);
                v_interp(i) += wk*Ncoeff[k][k];
            }                        

        }
        if(i>0)
            last_ds = s1-s0;     
        else
            last_ds = tf(1)-tf(0);
    }


    return v_interp;
}

Vector Boundary::interpolate_nonsmooth(const Vector& v, const Boundary& bnd) const
{    
    Vector v_interp = Vector::Zero(n_pts);
    Real tol = 1e-13;
    
    const Vector& tc = bnd.t;
    const Vector& tf = this->t;
    
    unsigned n_pts_c = bnd.n_pts;
    unsigned n_pts_f = this->n_pts;
    
    Vector dtc(n_pts_c);
    dtc.segment(0,n_pts_c-1) = tc.segment(1,n_pts_c-1)-tc.segment(0,n_pts_c-1);
    dtc(n_pts_c-1) = 1.-tc(n_pts_c-1)+tc(0);
    
    Vector dtf(n_pts_f);
    dtf.segment(0,n_pts_f-1) = tf.segment(1,n_pts_f-1)-tf.segment(0,n_pts_f-1);
    dtf(n_pts_f-1) = 1.-tf(n_pts_f-1)+tf(0);        
    
    unsigned j1 = n_pts_c-1;
    Real t1 = tc(0)-dtc(n_pts_c-1);
    for(unsigned i=0;i<n_pts_f;i++)
    {
        while(t1+dtc(j1)<=tf(i)+tol/2.)
        {
            t1 += dtc(j1);
            j1 = (j1+1)%n_pts_c;                      
        }
        
        if(abs(tf(i)-t1)<tol)
            v_interp(i)=v(j1);
        else
        {            
            unsigned j2=(j1+1)%n_pts_c;
            Real t2 = t1 + dtc(j1);                        
            
            unsigned j3;
            Real t3;
            if(abs(tf(i)-t1)<abs(tf(i)-t2))
            {
                j3 = (j1-1+n_pts_c)%n_pts_c;
                t3 = t1-dtc(j3); 
            }
            else
            {
                j3 = (j2+1)%n_pts_c;
                t3 = t2 + dtc(j2);
            }
            Real dt1 = t1-tf(i);
            Real dt2 = t2-tf(i);
            Real dt3 = t3-tf(i);
            Real a1 = (dt3-dt2)/dt1;
            Real a2 = (dt1-dt3)/dt2;
            Real a3 = (dt2-dt1)/dt3;
            v_interp(i) = (a1*v(j1)+a2*v(j2)+a3*v(j3))/(a1+a2+a3);
//            Real a1 = -dt2/(dt1-dt2);
//            Real a2 = dt1/(dt1-dt2);
//            v_interp(i) = (a1*v(j1)+a2*v(j2))/(a1+a2);
//            cout<<"sum a = "<<a1+a2<<endl;
        }
        
    }
    
    return v_interp;
}

Real Boundary::integrate(const Vector& v) const
{
    if(v.size()!=n_pts)
    {
        cerr<<"ERROR: in Boundary::integrate"<<endl;
        return 0.;
    }
        
    return (v.array()*measure.array()).sum();
}

Real Boundary::L2norm(const Vector& v) const
{
    if(v.size()!=n_pts)
    {
        cerr<<"ERROR: in Boundary::integrate"<<endl;
        return 0.;
    }
        
    return sqrt((v.array().pow(2)*measure.array()).sum());
}

void Boundary::save(string name) const
{
    ofstream ofile(name+".bin", std::ios::out | std::ios::binary | std::ios::trunc);
    save(ofile);
    ofile.close();
}

void Boundary::load(string name)
{
    std::ifstream ifile(name+".bin", std::ios::in | std::ios::binary);
    load(ifile);
    ifile.close();
}

void Boundary::save(std::ofstream& out) const
{
    write_scalar_binary(out,n_pts);
    write_scalar_binary(out,h);
    write_scalar_binary(out,length);
    
    write_binary(out,t);
    write_binary(out,dt);
    write_binary(out,G);
    write_binary(out,dG);
    write_binary(out,ddG);
    write_binary(out,N);
    write_binary(out,norm_dG);
    write_binary(out,norm_ddG);
    
    write_scalar_binary(out,smooth_boundary);
    write_binary(out,FGx.get_FCoeff());
    write_binary(out,FGy.get_FCoeff());
    
    write_binary<unsigned,3>(out,sides_dofs);
    write_binary<unsigned,2>(out,corners);
}

void Boundary::load(std::ifstream& in)
{   
    Vector FGxcoeff,FGycoeff;
    
    read_scalar_binary(in,n_pts);
    read_scalar_binary(in,h);
    read_scalar_binary(in,length);
    
    read_binary(in,t);
    read_binary(in,dt);
    read_binary(in,G);
    read_binary(in,dG);
    read_binary(in,ddG);
    read_binary(in,N);
    read_binary(in,norm_dG);
    read_binary(in,norm_ddG);
    
    read_scalar_binary(in,smooth_boundary);
    read_binary(in,FGxcoeff);
    read_binary(in,FGycoeff);
    
    read_binary<unsigned,3>(in,sides_dofs);
    read_binary<unsigned,2>(in,corners);
    
    FGx.set_FCoeff(FGxcoeff);
    FGy.set_FCoeff(FGycoeff);
    
    dFGx.set_FCoeff(FGxcoeff);
    dFGy.set_FCoeff(FGycoeff);
    dFGx.diff();
    dFGy.diff();
    
    ddFGx.set_FCoeff(FGxcoeff);
    ddFGy.set_FCoeff(FGycoeff);
    ddFGx.diff();
    ddFGy.diff();
    ddFGx.diff();
    ddFGy.diff();
    
    measure = norm_dG.array()*dt.array();
}

bool Boundary::operator==(const Boundary& bnd) const
{
    if(n_pts!=bnd.get_n_pts())
        return false;    
    if(h!=bnd.get_h())
        return false;
    if(length!=bnd.get_length())
        return false;
    
    if(t!=bnd.get_t())
        return false;
    if(dt!=bnd.get_dt())
        return false;
    if(G!=bnd.get_G())
        return false;
    if(dG!=bnd.get_dG())
        return false;
    if(ddG!=bnd.get_ddG())
        return false;
    if(N!=bnd.get_N())
        return false;
    if(norm_dG!=bnd.get_norm_dG())
        return false;
    if(norm_ddG!=bnd.get_norm_ddG())
        return false;    
    if(measure!=bnd.measure)
        return false;
    
    if(sides_dofs!=bnd.get_sides_dofs())
        return false;
    if(corners!=bnd.get_corners())
        return false;
    
    if(smooth_boundary!=bnd.smooth_boundary)
        return false;
    
    if(FGx.get_FCoeff()!=bnd.FGx.get_FCoeff())
        return false;
    if(FGy.get_FCoeff()!=bnd.FGy.get_FCoeff())
        return false;
    if(dFGx.get_FCoeff()!=bnd.dFGx.get_FCoeff())
        return false;
    if(dFGy.get_FCoeff()!=bnd.dFGy.get_FCoeff())
        return false;
    if(ddFGx.get_FCoeff()!=bnd.ddFGx.get_FCoeff())
        return false;
    if(ddFGy.get_FCoeff()!=bnd.ddFGy.get_FCoeff())
        return false;
    
    return true;
}

bool Boundary::operator!=(const Boundary& bnd) const
{
    return !operator==(bnd);
}

bool Boundary::operator%(const Boundary& bnd) const
{
    if(n_pts!=bnd.get_n_pts())
        return false;   

    Real tol=1e-8; 
    
    if((t-bnd.get_t()).norm()>tol*min(t.norm(),bnd.get_t().norm()))
        return false;
    
    // Verify that the two boundary are the same up to a translation
    Point p = G.row(0)-bnd.G.row(0);
    Real norm_p = p.norm();
    for(unsigned i=0;i<G.rows();i++)
        if((G.row(i)-bnd.G.row(i)-p).norm()>tol*norm_p)
            return false;

    //don't check the remaining...
    // if((dG-bnd.get_dG()).norm()>tol*min(dG.norm(),bnd.get_dG().norm()))
    //     return false;
    // if((ddG-bnd.get_ddG()).norm()>tol*min(ddG.norm(),bnd.get_ddG().norm()))
    //     return false;

    // if(N!=bnd.get_N())
    //     return false;
    // if(norm_dG!=bnd.get_norm_dG())
    //     return false;
    // if(norm_ddG!=bnd.get_norm_ddG())
    //     return false;    
    // if(measure!=bnd.measure)
    //     return false;
    
    // if(sides_dofs!=bnd.get_sides_dofs())
    //     return false;
    // if(corners!=bnd.get_corners())
    //     return false;
    
    // if(smooth_boundary!=bnd.smooth_boundary)
    //     return false;
    
    return true;
}

ostream& operator<<(ostream& os, const uBidomain::Boundary& bnd)
{
    auto n = bnd.get_n_pts();
    os<<"n = "<<n<<";"<<endl;
    
    write_matrix(os,bnd.get_t(),"t");
    write_matrix(os,bnd.get_G(),"G");
    write_matrix(os,bnd.get_dG(),"dG");
    write_matrix(os,bnd.get_ddG(),"ddG");
    write_matrix(os,bnd.get_N(),"N");
    write_matrix(os,bnd.get_norm_dG(),"norm_dG");
    write_matrix(os,bnd.get_norm_ddG(),"norm_ddG");   
    
    return os;
}
