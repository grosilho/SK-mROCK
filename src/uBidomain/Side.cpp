#include <vector>

#include "Side.h"
using namespace uBidomain;

Side::Side()
:x(0), G(0,0), dG(0,0), ddG(0,0), N(0,0), 
 norm_dG(0), norm_ddG(0),
 p1(0.,0.), p2(0.,0.), 
 n(0), dx(0.), length(0.), h(0.), kappa(0)
{
}

Side::Side(const vector<Side>& sides)
{
    concatenate_and_set(sides);
}

void Side::concatenate_and_set(const vector<Side>& sides)
{
    //Set this side as concatenation of the given ones
    //the given sides must be ordered
    
    unsigned n_sides = sides.size();
    
    n=0;
    for(auto side:sides)
        n += side.get_n();
    
    x = Eigen::ArrayXd::LinSpaced(n, 0, n - 1) / n;
    dx = 1./n;
    x = x.array()+dx/2.;
    
    p1 = sides.front().begin();
    p2 = sides.back().end();
    
    G = Matrix::Zero(n,2);
    dG = Matrix::Zero(n,2);
    ddG = Matrix::Zero(n,2);
    N = Matrix::Zero(n,2);
    norm_dG = Vector::Zero(n);
    norm_ddG = Vector::Zero(n);
    h = std::numeric_limits<Real>::max();
    length = 0.;
    
    kappa = Vector::Zero(n);
        
    unsigned dofs_count=0;
    
    Real tol=1e-10;
    
    for(auto i=0u;i<n_sides;i++)
    {        
        const Side& side = sides[i];        
        unsigned side_size = side.get_n();
        
        Real dx_side = side.get_dx();
        G.block(dofs_count,0,side_size,2) = side.get_G();
        dG.block(dofs_count,0,side_size,2) = (dx_side/dx)*side.get_dG();
        ddG.block(dofs_count,0,side_size,2) = (dx_side/dx)*(dx_side/dx)*side.get_ddG();
        N.block(dofs_count,0,side_size,2) = side.get_N();
        norm_dG.segment(dofs_count,side_size) = (dx_side/dx)*side.get_norm_dG();
        norm_ddG.segment(dofs_count,side_size) = (dx_side/dx)*(dx_side/dx)*side.get_norm_ddG();
        h = min(h,side.get_h());
        length += side.get_length();
                        
        kappa.segment(dofs_count,side_size) = side.get_kappa();
        
        dofs_count += side_size;
        
        if(i<n_sides-1)
            if((sides[i].end()-sides[i+1].begin()).norm()>tol)
                cout<<"ERROR in building concatenated side."<<endl;     
    }
}

void Side::resize_vectors(unsigned int n)
{
    x.resize(n);
    G.resize(n,2);
    dG.resize(n,2);
    ddG.resize(n,2);
    N.resize(n,2);
    norm_dG.resize(n);
    norm_ddG.resize(n);
    kappa.resize(n);
}

Side Side::reverse_normal()
{
    Side S(*this);
    S.N = -S.N;
    
    return S;
}

Side Side::operator-() const
{
    Side mS(*this);
    
    mS.x = this->x(Eigen::seq(Eigen::last,0,Eigen::fix<-1>));
    mS.G = this->G(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
    mS.dG = -this->dG(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
    mS.ddG = this->ddG(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
    mS.N = -this->N(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
    mS.norm_dG = this->norm_dG(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
    mS.norm_ddG = this->norm_ddG(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
    mS.p1 = this->p2;
    mS.p2 = this->p1;
    
    mS.kappa = this->kappa(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);    
    
    return mS;
}

Side Side::operator+(const Point& p) const
{
    Side S(*this);
    
    S.G.rowwise() += p;;
    S.p1 += p;
    S.p2 += p;
    
    return S;
}

Side Side::operator*(const Real r) const
{
    Side S(*this);
    
    S.G = S.G*r;
    S.dG = S.dG*r;
    S.ddG = S.ddG*r;
    S.norm_dG = S.norm_dG*r;
    S.norm_ddG = S.norm_ddG*r;
    S.p1 = S.p1*r;
    S.p2 = S.p2*r;
    S.length = S.length*r;
    S.h = S.h*r;
    
    return S;
}

const Vector& Side::get_x() const
{
    return x;
}

const Matrix& Side::get_G() const
{
    return G;
}

const Matrix& Side::get_dG() const
{
    return dG;
}

const Matrix& Side::get_ddG() const
{
    return ddG;
}

const Matrix& Side::get_N() const
{
    return N;
}

const Vector& Side::get_norm_dG() const
{
    return norm_dG;
}

const Vector& Side::get_norm_ddG() const
{
    return norm_ddG;
}

const Vector& Side::get_kappa() const
{
    return kappa;
}

Point Side::get_G(unsigned int i) const
{
    return G.block<1,2>(i,0);
}

Point Side::get_N(unsigned int i) const
{
    return N.block<1,2>(i,0);
}

Real Side::get_norm_dG(unsigned int i) const
{
    return norm_dG(i);
}

unsigned int Side::get_n() const
{
    return n;
}

Real Side::get_dx() const
{
    return dx;
}

Real Side::get_h() const
{
    return h;
}

Real Side::get_length() const
{
    return length;
}

Point Side::begin() const
{
    return p1;
}

Point Side::end() const
{
    return p2;
}

ostream& Side::write(ostream& os, string name)
{
    if(!name.empty())
        name += ".";
    
    os<<name<<"n = "<<n<<";"<<endl;
    
    write_matrix(os,x,name+"x");
    write_matrix(os,G,name+"G");
    write_matrix(os,dG,name+"dG");
    write_matrix(os,ddG,name+"ddG");
    write_matrix(os,N,name+"N");
    write_matrix(os,norm_dG,name+"norm_dG");
    write_matrix(os,norm_ddG,name+"norm_ddG");   
    write_matrix(os,kappa,name+"kappa");   
    
    return os;
}
// -----------------------------------------------------------------------------

ostream& operator<<(ostream& os, const Matrix& M)
{
    auto n = M.rows();
    auto m = M.cols();
    
    os<<"[";
    for(unsigned i=0;i<n-1;i++)
    {
        for(unsigned j=0;j<m-1;j++)
            os<<M(i,j)<<",";
        os<<M(i,m-1)<<";";
    }
    for(unsigned j=0;j<m-1;j++)
        os<<M(n-1,j)<<",";
    os<<M(n-1,m-1)<<"]"<<endl;
    
    return os;
}

void write_matrix(ostream& os, const Matrix& M, string name)
{
    auto n = M.rows();
    auto m = M.cols();
    
    if(n==0 || m==0)
        return;
    
    os<<name<<" = zeros("<<n<<","<<m<<");"<<endl;
    
    for(unsigned j=0;j<m;j++)
    {
        os<<name<<"(:,"<<j+1<<")=[";
        for(unsigned i=0;i<n-1;i++)
            os<<M(i,j)<<";";
        os<<M(n-1,j)<<"];"<<endl;
    }
}

ostream& operator<<(ostream& os, const Side& s)
{
    auto n = s.get_n();
    os<<"n = "<<n<<";"<<endl;
    
    write_matrix(os,s.get_x(),"x");
    write_matrix(os,s.get_G(),"G");
    write_matrix(os,s.get_dG(),"dG");
    write_matrix(os,s.get_ddG(),"ddG");
    write_matrix(os,s.get_N(),"N");
    write_matrix(os,s.get_norm_dG(),"norm_dG");
    write_matrix(os,s.get_norm_ddG(),"norm_ddG");   
    write_matrix(os,s.get_kappa(),"kappa");   
    
    return os;
}

Side operator*(const Real r,const Side& s)
{
    return s*r;
}

Side operator+(const Point p,const Side& s)
{
    return s+p;
}

Side operator-(const Point p,const Side& s)
{
    return (-s)+p;
}

// -----------------------------------------------------------------------------

Segment::Segment(Point p1_, Point p2_, Real h_, Real sign_N, Real kappa_)
:Side()
{   
    p1 = p1_;
    p2 = p2_;
    length = (p2-p1).norm();
    n = ceil(length/h_);
    n = n + n%2; //make sure it is even
    dx = 1./n;
    h = length/n;
    
    resize_vectors(n);
    
    kappa.array() = kappa_;
    
    Point tangent((p2-p1)/length);
    Point normal(tangent(1),-tangent(0));
    normal *= sign_N;
    for(unsigned int i=0;i<n;i++)
    {   
        x(i) = (2.*i+1.)/(2.*n);
        G.block(i,0,1,2) = (p1 + (p2-p1)*(2*i+1.)/(2.*n));
        dG.block(i,0,1,2) = p2-p1;
        ddG.block(i,0,1,2) = Point{0.,0.};
        N.block(i,0,1,2) = normal;
        norm_dG(i) = length;
        norm_ddG(i) = 0.;
    }
}



// -----------------------------------------------------------------------------

Arc::Arc(Point c_, Real r_, Real alpha_, Real beta_, Real h_, Real sign_N, Real kappa_)
:Side(), c(c_), r(r_), alpha(alpha_), beta(beta_)
{   
    p1 = c + r*Point(cos(alpha),sin(alpha));
    p2 = c + r*Point(cos(beta),sin(beta));
    
    length = abs(beta-alpha)*r;
    n = ceil(length/h_);
    n = n + n%2;
    dx = 1./n;
    h = length/n;
    
    resize_vectors(n);
    
    kappa.array() = kappa_;
    
    Real angle;
    for(unsigned int i=0;i<n;i++)
    {   
        x(i) = (2.*i+1.)/(2.*n);
        angle = alpha+(beta-alpha)*x(i);
        G.block(i,0,1,2) = c + r*Point(cos(angle),sin(angle));
        dG.block(i,0,1,2) = r*(beta-alpha)*Point(-sin(angle),cos(angle));
        ddG.block(i,0,1,2) = r*(beta-alpha)*(beta-alpha)*Point(-cos(angle),-sin(angle));
        N.block(i,0,1,2) = sign_N*Point(cos(angle),sin(angle));
        norm_dG(i) = r*abs(beta-alpha);
        norm_ddG(i) = r*(beta-alpha)*(beta-alpha);
    }
}

// -----------------------------------------------------------------------------

EllipseArc::EllipseArc(Point c_, Real r1_, Real r2_, 
                       Real alpha_, Real beta_, Real h_, Real sign_N, Real kappa_)
:Side(), c(c_), r1(r1_), r2(r2_), alpha(alpha_), beta(beta_)
{   
    p1 = c + Point(r1*cos(alpha),r2*sin(alpha));
    p2 = c + Point(r1*cos(beta),r2*sin(beta));
    
    //compute perimeter numerically, there's no formula for ellipse
    auto norm_dg = [=](Real t){return sqrt(pow(r1*sin(t),2)+pow(r2*cos(t),2));};
    length =0.;
    Real dtheta=0.001;
    unsigned Ntheta = round(abs(beta-alpha)/dtheta);
    dtheta = (beta-alpha)/Ntheta;
    for(unsigned i=0;i<Ntheta;i++)
        length += norm_dg(alpha+i*dtheta)*abs(dtheta);
    
    n = ceil(length/h_);
    n = n + n%2;
    dx = 1./n;
    h = length/n;
    
    resize_vectors(n);
    
    kappa.array() = kappa_;
    
    Real angle;
    for(unsigned int i=0;i<n;i++)
    {   
        x(i) = (2.*i+1.)/(2.*n);
        angle = alpha+(beta-alpha)*x(i);
        G.block(i,0,1,2) = c + Point(r1*cos(angle),r2*sin(angle));
        dG.block(i,0,1,2) = (beta-alpha)*Point(-r1*sin(angle),r2*cos(angle));
        ddG.block(i,0,1,2) = (beta-alpha)*(beta-alpha)*Point(-r1*cos(angle),-r2*sin(angle));
        N.block(i,0,1,2) = sign_N*Point(r2*cos(angle),r1*sin(angle)).normalized();
        norm_dG(i) = abs(beta-alpha)*Point(-r1*sin(angle),r2*cos(angle)).norm();
        norm_ddG(i) = (beta-alpha)*(beta-alpha)*Point(-r1*cos(angle),-r2*sin(angle)).norm();
    }
}

// -----------------------------------------------------------------------------

Wave1::Wave1(Point p1_, Point p2_, unsigned k_, Real a_, Real h_, Real sign_N, Real kappa_)
:Side(), k(k_), a(a_)
{   
    p1=p1_;
    p2=p2_;
    
    //compute length numerically
    Real p1p2sq = (p2-p1).squaredNorm();
    auto norm_dg = [=](Real t){return sqrt( p1p2sq+pow(2.*a*k*M_PI*sin(t*k*2.*M_PI),2) );};
    length =0.;
    Real dt=min(0.0001,0.01/k);
    unsigned Ndt = round(1./dt);
    dt = 1./Ndt;
    for(unsigned i=0;i<Ndt;i++)
        length += norm_dg((i+0.5)*dt);
    length*=dt;
    
    //minimum mesh size h needed to capture the wave
    if(k>0.)
        h_ = min(h_,length/k/32.);
    
    n = ceil(length/h_);
    n = n + n%2;
    dx = 1./n;
    h = length/n;
        
    resize_vectors(n);
    
    kappa.array() = kappa_;
    
    Point p2p1 = p2-p1;
    Point N_p2p1 = Point(p2p1(1),-p2p1(0));//rotated pi/2 clockwise
    N_p2p1 = N_p2p1/N_p2p1.norm();
    
    for(unsigned int i=0;i<n;i++)
    {   
        x(i) = (2.*i+1.)/(2.*n);
        G.block(i,0,1,2) = p1+p2p1*x(i)+a*(1.-cos(x(i)*2.*k*M_PI))*N_p2p1;
        dG.block(i,0,1,2) = p2p1+2.*a*k*M_PI*sin(x(i)*2.*k*M_PI)*N_p2p1;
        ddG.block(i,0,1,2) = -4.*k*M_PI*a*k*M_PI*cos(x(i)*2.*k*M_PI)*N_p2p1;
        N.block(i,0,1,2) = sign_N*Point(dG(i,1),-dG(i,0)).normalized();
        norm_dG(i) = norm_dg(x(i));
        norm_ddG(i) = 4.*k*M_PI*a*k*M_PI*abs(cos(x(i)*2.*k*M_PI));
    }
}

// -----------------------------------------------------------------------------

Wave2::Wave2(Point p1_, Point p2_, unsigned k_, Real a_, Real h_, Real sign_N, Real kappa_)
:Side(), k(k_), a(a_)
{   
    p1=p1_;
    p2=p2_;
    
    //compute length numerically
    Real p1p2sq = (p2-p1).squaredNorm();
    auto norm_dg = [=](Real t){return sqrt( p1p2sq+pow(2.*a*k*M_PI*cos(t*k*2.*M_PI),2) );};
    length =0.;
    Real dt=min(0.0001,0.01/k);
    unsigned Ndt = round(1./dt);
    dt = 1./Ndt;
    for(unsigned i=0;i<Ndt;i++)
        length += norm_dg((i+0.5)*dt);
    length*=dt;
    
    //minimum mesh size h needed to capture the wave
    if(k>0.)
        h_ = min(h_,length/k/32.);
    
    n = ceil(length/h_);
    n = n + n%2;
    dx = 1./n;
    h = length/n;
        
    resize_vectors(n);
    
    kappa.array() = kappa_;
    
    Point p2p1 = p2-p1;
    Point N_p2p1 = Point(p2p1(1),-p2p1(0));//rotated pi/2 clockwise
    N_p2p1 = N_p2p1/N_p2p1.norm();
    
    for(unsigned int i=0;i<n;i++)
    {   
        x(i) = (2.*i+1.)/(2.*n);
        G.block(i,0,1,2) = p1+p2p1*x(i)+a*sin(x(i)*2.*k*M_PI)*N_p2p1;
        dG.block(i,0,1,2) = p2p1+2.*a*k*M_PI*cos(x(i)*2.*k*M_PI)*N_p2p1;
        ddG.block(i,0,1,2) = -4.*k*M_PI*a*k*M_PI*sin(x(i)*2.*k*M_PI)*N_p2p1;
        N.block(i,0,1,2) = sign_N*Point(dG(i,1),-dG(i,0)).normalized();
        norm_dG(i) = norm_dg(x(i));
        norm_ddG(i) = 4.*k*M_PI*a*k*M_PI*abs(sin(x(i)*2.*k*M_PI));
    }
}

// -----------------------------------------------------------------------------

SquaredWave1::SquaredWave1(Point p1_, Point p2_, unsigned k_, Real a_, Real h_, Real sign_N, Real kappa_dir, Real kappa_perp)
:Side(), k(k_), a(a_)
{   
    p1=p1_;
    p2=p2_;
    h=h_;
    
//    Point dir = p2-p1;
//    Real ls = dir.norm()/k/4.;
//    Real hs = a/2.;
//    dir = dir/dir.norm();
//    Point perp(dir(1),-dir(0));
//    
//    Side h_half = Segment({0.,0.},(ls/2.)*dir,h);
//    Side h_full = Segment({0.,0.},ls*dir,h);
//    Side v = Segment({0.,0.},hs*perp,h);
//    Side vr = Segment({0.,0.},-hs*perp,h);
//    
//    vector<Side> sides(8*k+1);
//    sides[0] = p1+h_half;
//    for(unsigned j=0;j<k;j++)
//    {     
//        sides[j*8+1] = sides[j*8+0].end()+v;
//        sides[j*8+2] = sides[j*8+1].end()+h_full;
//        sides[j*8+3] = sides[j*8+2].end()+v;
//        sides[j*8+4] = sides[j*8+3].end()+h_full;
//        sides[j*8+5] = sides[j*8+4].end()+vr;
//        sides[j*8+6] = sides[j*8+5].end()+h_full;
//        sides[j*8+7] = sides[j*8+6].end()+vr;
//        if(j<k-1)
//            sides[j*8+8] = sides[j*8+7].end()+h_full;
//    }
//    sides[k*8] = sides[k*8-1].end()+h_half;
//    
//    concatenate_and_set(sides);

    Point dir = p2-p1;
    Real ls = dir.norm()/k/2.;
    Real hs = a;
    dir = dir/dir.norm();
    Point perp(dir(1),-dir(0));
    
    Side h_half = Segment({0.,0.},(ls/2.)*dir,h,1.,kappa_dir);
    Side h_full = Segment({0.,0.},ls*dir,h,1.,kappa_dir);
    Side v = Segment({0.,0.},hs*perp,h,1.,kappa_perp);
    Side vr = Segment({0.,0.},-hs*perp,h,1.,kappa_perp);
    
    vector<Side> sides(4*k+1);
    sides[0] = p1+h_half;
    for(unsigned j=0;j<k;j++)
    {     
        sides[j*4+1] = sides[j*4+0].end()+v;
        sides[j*4+2] = sides[j*4+1].end()+h_full;
        sides[j*4+3] = sides[j*4+2].end()+vr;        
        if(j<k-1)
            sides[j*4+4] = sides[j*4+3].end()+h_full;
    }
    sides[k*4] = sides[k*4-1].end()+h_half;
    
    concatenate_and_set(sides);
}

// -----------------------------------------------------------------------------

SquaredWave2::SquaredWave2(Point p1_, Point p2_, unsigned k_, Real a_, Real h_, Real sign_N, Real kappa_dir, Real kappa_perp)
:Side(), k(k_), a(a_)
{   
    p1=p1_;
    p2=p2_;
    h=h_;
    
    Point dir = p2-p1;
    Real ls = dir.norm()/k/2.;
    Real hs = 2.*a;
    dir = dir/dir.norm();
    Point perp(dir(1),-dir(0));
    
    Side h_full = Segment({0.,0.},ls*dir,h,1.,kappa_dir);
    Side v = Segment({0.,0.},hs*perp,h,1.,kappa_perp);
    Side v_half = Segment({0.,0.},0.5*hs*perp,h,1.,kappa_perp);
    Side vr = Segment({0.,0.},-hs*perp,h,1.,kappa_perp);
    
    vector<Side> sides(4*k+1);
    sides[0] = p1+v_half;
    for(unsigned j=0;j<k;j++)
    {     
        sides[j*4+1] = sides[j*4+0].end()+h_full;
        sides[j*4+2] = sides[j*4+1].end()+vr;
        sides[j*4+3] = sides[j*4+2].end()+h_full;          
        if(j<k-1)
            sides[j*4+4] = sides[j*4+3].end()+v;
    }
    sides[k*4] = sides[k*4-1].end()+v_half;
    
    concatenate_and_set(sides);
}