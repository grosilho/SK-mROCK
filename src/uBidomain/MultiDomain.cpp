#include <fstream>
#include "MultiDomain.h"
#include "Utils.h"
#include <random>
#include <chrono>

using namespace uBidomain;

void add_cells_cluster(const Parameters& param, Point& p, 
        vector<SidesToBoundaryMap>& cells, vector<Side>& sides);
//void add_cells_cluster_with_holes(const Real cl, const Real cw, 
//        const Real pl, const Real pw, Point& p, 
//        unsigned nx, unsigned ny, const Real dGl, const Real dGw,
//        vector<SidesToBoundaryMap>& cells, vector<Side>& sides);
//Side squared_piramid(const Point a, const Point b, 
//                     const Real H, const Real dG);

MultiDomain::MultiDomain()
{
    
}

MultiDomain::MultiDomain(const Boundary& outer_bnd_, const vector<Side>& sides_, 
                const vector<SidesToBoundaryMap> inner_dom_desc)
{
    init_MultiDomain(outer_bnd_,sides_,inner_dom_desc);
}

void MultiDomain::init_MultiDomain(const Boundary& outer_bnd_, const vector<Side>& sides_, 
                const vector<SidesToBoundaryMap> inner_dom_desc)
{    
    outer_bnd=outer_bnd_;
    sides=sides_;
    
    // inner domains plus the outer domain
    cout<<"Building MultiDomain..."<<endl;
    
    N_domains = 1+inner_dom_desc.size(); 
    N_inner_domains = inner_dom_desc.size(); 
    N_sides = sides.size();
    N_outer_bnd_dofs = outer_bnd.get_n_pts();
    
    // we define the global index and size of each side
    global_dofs.resize(N_sides);
    unsigned dofs_count=0;
    for(auto i=0u;i<sides.size();i++)
    {
        global_dofs[i] = {dofs_count, sides[i].get_n()};
        dofs_count += sides[i].get_n();
    }
    N_global_dofs = dofs_count;
    
    //This vector counts how many times a side is used to construct the inner 
    //boundaries. If a side is used only once then it is part of the outer domain.
    //If it is used twice then it is a gap junction.
    vector<unsigned> count_side(sides.size(),0);
    //For every side we get the index of its neighboring domains
    vector<vector<int>> neigh_doms(sides.size());
   
    //for every domain we get a list of
    //-each side in the domain, represented by its index. The index is the same used to access global_dofs
    //-for each side we get the direction, its local dofs and its size (the size is the same as the one in global_dofs).
    doms_GtL_maps.resize(N_domains);
    //we keep track of the number of dofs in each domain
    N_local_dofs.resize(N_domains);
    
    //We reserve space for the outer domain. Will be corrected later.
    doms.clear();
    doms.push_back(Domain(outer_bnd)); 
    
    // build the inner domains
    for(unsigned i=0;i<inner_dom_desc.size();i++)
    {
        // take the object describing this inner domain
        // it contains a list of sides (their indeces) and their direction
        const SidesToBoundaryMap& StDm = inner_dom_desc[i];
        
        unsigned n_sides = StDm.size();
        vector<Side> dom_sides(n_sides);
        for(unsigned j=0;j<n_sides;j++)
        {
            if(StDm[j][1]==1)// use the side in the direction it is
                dom_sides[j] = sides[StDm[j][0]];
            else//use the side in the reverse direction (reverse the nodes order and change the normal sign)
                dom_sides[j] = -sides[StDm[j][0]];
            
            //keep track of how many times we use this side
            count_side[StDm[j][0]]++;
            //and which domain has this side as part of its boundary
            neigh_doms[StDm[j][0]].push_back(i+1); //+1 because the first domain is the outer domain
        }
        //build the domain. Note that the list of sides dom_sides can be unordered.
        //The constructor of Boundary will sort them in counter clock wise order.
        //However, the normal must be already in the good direction.
        doms.push_back(Domain(Boundary(dom_sides)));
        
        //now we define the sides local degrees of freedom
        GlobalToLocalMap this_dom_loc_dofs(n_sides);
        //for each side, sides_dofs contains the starting index and length
        //of this side. 
        //side_to_dofs is ordered in the counter clock wise sense defined by
        //the constructor of Boundary. However it contains info to recover the previous order
        vector<array<unsigned,3>> sides_dofs 
            = doms.back().get_bnds()[0].get_sides_dofs();
        for(unsigned j=0;j<n_sides;j++)
        {
            unsigned prev_j = sides_dofs[j][2]; //original order given to constructor of Boundary
            this_dom_loc_dofs[j][0] = StDm[prev_j][0]; //the index of the side in the global dof system. 
            this_dom_loc_dofs[j][1] = global_dofs[StDm[prev_j][0]][0]; //starting index of the side in the global dof system                    
            this_dom_loc_dofs[j][2] = sides_dofs[j][0]; //starting index of the side in the local dof system
            this_dom_loc_dofs[j][3] = sides_dofs[j][1]; //length of the side
            this_dom_loc_dofs[j][4] = StDm[prev_j][1]; //direction of the side
            // it remains to define this_dom_loc_dofs[j][5], which tells if in B_i we must do a sign change or not. It's done later.
        }
        doms_GtL_maps[i+1] = this_dom_loc_dofs; //+1 because first one is outer domain
        N_local_dofs[i+1] = doms.back().get_tot_n_dofs();
    }
    
    //now we take the sides used only once and define the inner boundary of the outer domain
    vector<Side> inner_sides;
    vector<unsigned> inner_sides_ind;
    for(auto i=0u;i<sides.size();i++)
        if(count_side[i]==1)
        {
            inner_sides.push_back(sides[i].reverse_normal());
            inner_sides_ind.push_back(i);
            neigh_doms[i].push_back(0);// side i touches domain 0
        }
    
    Real tol=1e-10;
    vector<vector<Side>> boundaries;
    vector<vector<unsigned>> sides_indeces;
    boundaries.clear();
    sides_indeces.clear();
    bool found;
    
    while(!inner_sides.empty())
    {
        boundaries.push_back(vector<Side>(0));
        boundaries.back().push_back(inner_sides.front());
        inner_sides.erase(inner_sides.begin());
        
        sides_indeces.push_back(vector<unsigned>(0));
        sides_indeces.back().push_back(inner_sides_ind.front());
        inner_sides_ind.erase(inner_sides_ind.begin());
        
        do
        {
            found=false;
            for(unsigned i=0;i<inner_sides.size();i++)
            {
                if((boundaries.back().back().end()-inner_sides[i].begin()).norm()<tol)
                {
                    found=true;
                    boundaries.back().push_back(inner_sides[i]);
                    inner_sides.erase(inner_sides.begin()+i);
                    sides_indeces.back().push_back(inner_sides_ind[i]);
                    inner_sides_ind.erase(inner_sides_ind.begin()+i);
                    break;
                }
            }
        }while(found);
    }
    
    vector<Boundary> inner_bnds;
    for(auto& boundary_sides:boundaries)
        inner_bnds.push_back(Boundary(boundary_sides));
    
    //and get its local dofs
    GlobalToLocalMap this_dom_loc_dofs;
    unsigned prev_bnd_pts=0;
    for(unsigned i=0;i<inner_bnds.size();i++)
    {
        vector<array<unsigned,3>> sides_dofs = inner_bnds[i].get_sides_dofs();        
        for(unsigned j=0;j<sides_dofs.size();j++)
        {            
            unsigned prev_j = sides_dofs[j][2];
            this_dom_loc_dofs.push_back(array<int,6>());
            this_dom_loc_dofs.back()[0] = sides_indeces[i][prev_j]; //the index of the side in the global dof system. 
            this_dom_loc_dofs.back()[1] = global_dofs[sides_indeces[i][prev_j]][0]; //starting index of the side in the global dof system        
            this_dom_loc_dofs.back()[2] = prev_bnd_pts+sides_dofs[j][0]; //starting index of the side in the local domain
            this_dom_loc_dofs.back()[3] = sides_dofs[j][1]; //length of the side
            this_dom_loc_dofs.back()[4] = 1; //direction of the side
            this_dom_loc_dofs.back()[5] = -1;//the index of the outer domain is always smaller than the index of any other domain, so there is sign change in B_i
        }
        prev_bnd_pts+=inner_bnds[i].get_n_pts();
    }
    doms_GtL_maps[0] = this_dom_loc_dofs;
    N_local_dofs[0] = prev_bnd_pts;
    
    //define the outer domain
    inner_bnds.insert(inner_bnds.begin(),outer_bnd);
    doms[0]=Domain(inner_bnds);
    
    //detect the gap junctions and build the map Global->Gap Junctions
    gap_junct_GtL_map.clear();
    unsigned local_ind=0;
    for(auto i=0u;i<sides.size();i++)
        if(count_side[i]==2)
        {           
            gap_junct_GtL_map.push_back({(int)i,(int)global_dofs[i][0],(int)local_ind,(int)global_dofs[i][1],1,1});
            local_ind += global_dofs[i][1];
        }
    if(local_ind>0)
        N_gap_junct_dofs = gap_junct_GtL_map.back()[2]+gap_junct_GtL_map.back()[3];
    else
        N_gap_junct_dofs=0;
    N_transmembrane_dofs = N_global_dofs-N_gap_junct_dofs;
    //It remains to tell to each GlobaltoLocalMap if in this domain the side has
    //positive or negative sign. Depending on the index of this domain and the
    //neighboring one.
    for(auto i=1u;i<N_domains;i++)
    {
        GlobalToLocalMap& GtLmap = doms_GtL_maps[i];
        for(auto& side:GtLmap)
        {
            int ind = side[0];
            if(i==max(neigh_doms[ind][0],neigh_doms[ind][1]))
                side[5]=1;
            else if(i==min(neigh_doms[ind][0],neigh_doms[ind][1]))
                side[5]=-1;
            else
                cerr<<"ERROR in MultiDomain constructor"<<endl;
        }
    }
    
    h_min = h_min = std::numeric_limits<Real>::max();
    h_max = 0.;
    for(auto dom:doms)
    {
        h_min = min(h_min,dom.get_h_min());
        h_max = max(h_max,dom.get_h_max());
    }
}

void MultiDomain::add_AiTViAi(Matrix& V, const Matrix& Vi, const Real alpha, unsigned i, bool sign_change) const
{
    // This method performs V += alpha*A_i^T*Vi*A_i efficiently
    // If sign_change==true then it uses B_i instead of A_i
    
    for(unsigned l=0;l<doms_GtL_maps[i].size();l++)
    {
        int sign_l = sign_change ? doms_GtL_maps[i][l][5] : 1;                
        
        for(unsigned c=0;c<doms_GtL_maps[i].size();c++)
        {
            int sign_c = sign_change ? doms_GtL_maps[i][c][5] : 1;            
            
            if(doms_GtL_maps[i][l][4]==1 && doms_GtL_maps[i][c][4]==1)
                V.block(doms_GtL_maps[i][l][1],doms_GtL_maps[i][c][1],
                        doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3]) 
                    += alpha*sign_l*sign_c
                      *Vi.block(doms_GtL_maps[i][l][2],doms_GtL_maps[i][c][2],
                                doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])
                        (Eigen::all, Eigen::all);
            else if(doms_GtL_maps[i][l][4]==0 && doms_GtL_maps[i][c][4]==1)
                V.block(doms_GtL_maps[i][l][1],doms_GtL_maps[i][c][1],
                        doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3]) 
                    += alpha*sign_l*sign_c
                      *Vi.block(doms_GtL_maps[i][l][2],doms_GtL_maps[i][c][2],
                                doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])
                        (Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
            else if(doms_GtL_maps[i][l][4]==1 && doms_GtL_maps[i][c][4]==0)
                V.block(doms_GtL_maps[i][l][1],doms_GtL_maps[i][c][1],
                        doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3]) 
                    += alpha*sign_l*sign_c
                      *Vi.block(doms_GtL_maps[i][l][2],doms_GtL_maps[i][c][2],
                                doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])
                        (Eigen::all,Eigen::seq(Eigen::last,0,Eigen::fix<-1>));
            else if(doms_GtL_maps[i][l][4]==0 && doms_GtL_maps[i][c][4]==0)
                V.block(doms_GtL_maps[i][l][1],doms_GtL_maps[i][c][1],
                        doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3]) 
                    += alpha*sign_l*sign_c
                      *Vi.block(doms_GtL_maps[i][l][2],doms_GtL_maps[i][c][2],
                                doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])
                        (Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::seq(Eigen::last,0,Eigen::fix<-1>));                        
        }
    }
}

void MultiDomain::add_AiVAiT(Matrix& Vi, const Matrix& V, const Real alpha, unsigned i, bool sign_change) const
{
    // This method performs Vi += alpha*A_i*V*A_i^T efficiently
    // If sign_change==true then it uses B_i instead of A_i
    
    for(unsigned l=0;l<doms_GtL_maps[i].size();l++)
    {
        int sign_l = sign_change ? doms_GtL_maps[i][l][5] : 1;                
        
        for(unsigned c=0;c<doms_GtL_maps[i].size();c++)
        {
            int sign_c = sign_change ? doms_GtL_maps[i][c][5] : 1;            
            
            if(doms_GtL_maps[i][l][4]==1 && doms_GtL_maps[i][c][4]==1)
                 Vi.block(doms_GtL_maps[i][l][2],doms_GtL_maps[i][c][2],
                          doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])                        
                    += alpha*sign_l*sign_c
                      *V.block(doms_GtL_maps[i][l][1],doms_GtL_maps[i][c][1],
                               doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])
                              (Eigen::all, Eigen::all);
            else if(doms_GtL_maps[i][l][4]==0 && doms_GtL_maps[i][c][4]==1)
                Vi.block(doms_GtL_maps[i][l][2],doms_GtL_maps[i][c][2],
                          doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])                        
                    += alpha*sign_l*sign_c
                      *V.block(doms_GtL_maps[i][l][1],doms_GtL_maps[i][c][1],
                               doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])
                        (Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
            else if(doms_GtL_maps[i][l][4]==1 && doms_GtL_maps[i][c][4]==0)
                Vi.block(doms_GtL_maps[i][l][2],doms_GtL_maps[i][c][2],
                          doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])                        
                    += alpha*sign_l*sign_c
                      *V.block(doms_GtL_maps[i][l][1],doms_GtL_maps[i][c][1],
                               doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])
                        (Eigen::all,Eigen::seq(Eigen::last,0,Eigen::fix<-1>));
            else if(doms_GtL_maps[i][l][4]==0 && doms_GtL_maps[i][c][4]==0)
                Vi.block(doms_GtL_maps[i][l][2],doms_GtL_maps[i][c][2],
                          doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])                        
                    += alpha*sign_l*sign_c
                      *V.block(doms_GtL_maps[i][l][1],doms_GtL_maps[i][c][1],
                               doms_GtL_maps[i][l][3],doms_GtL_maps[i][c][3])
                        (Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::seq(Eigen::last,0,Eigen::fix<-1>));                        
        }
    }
}

void MultiDomain::add_AiVAgT(Matrix& Vi, const Matrix& V, const Real alpha, unsigned i, bool sign_change) const
{
    // This method performs Vi += alpha*A_i*V*A_g^T efficiently
    // If sign_change==true then it uses B_i instead of A_i, A_g remains the same
    
    for(unsigned l=0;l<doms_GtL_maps[i].size();l++)
    {
        int sign_l = sign_change ? doms_GtL_maps[i][l][5] : 1;                
        
        for(unsigned c=0;c<gap_junct_GtL_map.size();c++)
        {            
            if(doms_GtL_maps[i][l][4]==1)
                 Vi.block(doms_GtL_maps[i][l][2],gap_junct_GtL_map[c][2],
                          doms_GtL_maps[i][l][3],gap_junct_GtL_map[c][3])                        
                    += alpha*sign_l
                      *V.block(doms_GtL_maps[i][l][1],gap_junct_GtL_map[c][1],
                               doms_GtL_maps[i][l][3],gap_junct_GtL_map[c][3])
                              (Eigen::all, Eigen::all);
            else
                Vi.block(doms_GtL_maps[i][l][2],gap_junct_GtL_map[c][2],
                          doms_GtL_maps[i][l][3],gap_junct_GtL_map[c][3])                        
                    += alpha*sign_l
                      *V.block(doms_GtL_maps[i][l][1],gap_junct_GtL_map[c][1],
                               doms_GtL_maps[i][l][3],gap_junct_GtL_map[c][3])
                        (Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);            
        }
    }
}

void MultiDomain::add_AgVAiT(Matrix& Vi, const Matrix& V, const Real alpha, unsigned i, bool sign_change) const
{
    // This method performs Vg += alpha*A_g*V*A_i^T efficiently
    // If sign_change==true then it uses B_i instead of A_i
    
    for(unsigned l=0;l<gap_junct_GtL_map.size();l++)
    {        
        for(unsigned c=0;c<doms_GtL_maps[i].size();c++)
        {
            int sign_c = sign_change ? doms_GtL_maps[i][c][5] : 1;            
            
            if(doms_GtL_maps[i][c][4]==1)
                 Vi.block(gap_junct_GtL_map[l][2],doms_GtL_maps[i][c][2],
                          gap_junct_GtL_map[l][3],doms_GtL_maps[i][c][3])                        
                    += alpha*sign_c
                      *V.block(gap_junct_GtL_map[l][1],doms_GtL_maps[i][c][1],
                               gap_junct_GtL_map[l][3],doms_GtL_maps[i][c][3])
                              (Eigen::all, Eigen::all);
            else
                Vi.block(gap_junct_GtL_map[l][2],doms_GtL_maps[i][c][2],
                          gap_junct_GtL_map[l][3],doms_GtL_maps[i][c][3])                        
                    += alpha*sign_c
                      *V.block(gap_junct_GtL_map[l][1],doms_GtL_maps[i][c][1],
                               gap_junct_GtL_map[l][3],doms_GtL_maps[i][c][3])
                        (Eigen::all,Eigen::seq(Eigen::last,0,Eigen::fix<-1>));            
        }
    }
}

void MultiDomain::add_AgVAgT(Matrix& Vi, const Matrix& V, const Real alpha) const
{
    // This method performs Vg += alpha*A_g*V*A_g^T efficiently
    
    for(unsigned l=0;l<gap_junct_GtL_map.size();l++)
    {        
        for(unsigned c=0;c<gap_junct_GtL_map.size();c++)
        {            
             Vi.block(gap_junct_GtL_map[l][2],gap_junct_GtL_map[c][2],
                      gap_junct_GtL_map[l][3],gap_junct_GtL_map[c][3])                        
                += alpha
                  *V.block(gap_junct_GtL_map[l][1],gap_junct_GtL_map[c][1],
                           gap_junct_GtL_map[l][3],gap_junct_GtL_map[c][3]);        
        }
    }
}

Matrix MultiDomain::global_to_local_left(const Matrix& V, unsigned i, bool sign_change) const
{    
    Matrix Vi = Matrix::Zero(N_local_dofs[i],V.cols());
    for(unsigned j=0;j<doms_GtL_maps[i].size();j++)
    {
        int sign = sign_change ? doms_GtL_maps[i][j][5] : 1;
        if(doms_GtL_maps[i][j][4]==1)
            Vi.block(doms_GtL_maps[i][j][2],0,doms_GtL_maps[i][j][3],V.cols())=
                    sign*V.block(doms_GtL_maps[i][j][1],0,doms_GtL_maps[i][j][3],V.cols());
        else
        {
            Vi.block(doms_GtL_maps[i][j][2],0,doms_GtL_maps[i][j][3],V.cols())=
                    sign*V.block(doms_GtL_maps[i][j][1],0,doms_GtL_maps[i][j][3],V.cols())(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
        }
    }
    return Vi;
}

Matrix MultiDomain::global_to_local_right(const Matrix& Vi, unsigned i, bool sign_change) const
{
    Matrix V=Matrix::Zero(Vi.rows(),N_global_dofs);
    for(unsigned j=0;j<doms_GtL_maps[i].size();j++)
    {
        int sign = sign_change ? doms_GtL_maps[i][j][5] : 1;
        if(doms_GtL_maps[i][j][4]==1)
            V.block(0,doms_GtL_maps[i][j][1],Vi.rows(),doms_GtL_maps[i][j][3])=
                    sign*Vi.block(0,doms_GtL_maps[i][j][2],Vi.rows(),doms_GtL_maps[i][j][3]);
        else
        {
            V.block(0,doms_GtL_maps[i][j][1],Vi.rows(),doms_GtL_maps[i][j][3])=
                    sign*Vi.block(0,doms_GtL_maps[i][j][2],Vi.rows(),doms_GtL_maps[i][j][3])(Eigen::all,Eigen::seq(Eigen::last,0,Eigen::fix<-1>));
        }
    }
    return V;
}

Matrix MultiDomain::local_to_global_left(const Matrix& Vi, unsigned i, bool sign_change) const
{
    Matrix V = Matrix::Zero(N_global_dofs,Vi.cols());
    for(unsigned j=0;j<doms_GtL_maps[i].size();j++)
    {
        int sign = sign_change ? doms_GtL_maps[i][j][5] : 1;        
        if(doms_GtL_maps[i][j][4]==1)
            V.block(doms_GtL_maps[i][j][1],0,doms_GtL_maps[i][j][3],V.cols())=
                    sign*Vi.block(doms_GtL_maps[i][j][2],0,doms_GtL_maps[i][j][3],V.cols());
        else
        {
            V.block(doms_GtL_maps[i][j][1],0,doms_GtL_maps[i][j][3],V.cols())=
                    sign*Vi.block(doms_GtL_maps[i][j][2],0,doms_GtL_maps[i][j][3],V.cols())(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
        }
    }
    return V;
}

Matrix MultiDomain::local_to_global_right(const Matrix& V, unsigned i, bool sign_change) const
{
    Matrix Vi = Matrix::Zero(V.rows(),N_local_dofs[i]);
    for(unsigned j=0;j<doms_GtL_maps[i].size();j++)
    {
        int sign = sign_change ? doms_GtL_maps[i][j][5] : 1;
        if(doms_GtL_maps[i][j][4]==1)
            Vi.block(0,doms_GtL_maps[i][j][2],V.rows(),doms_GtL_maps[i][j][3])=
                    sign*V.block(0,doms_GtL_maps[i][j][1],V.rows(),doms_GtL_maps[i][j][3]);                  
        else
        {
            Vi.block(0,doms_GtL_maps[i][j][2],V.rows(),doms_GtL_maps[i][j][3])=
                    sign*V.block(0,doms_GtL_maps[i][j][1],V.rows(),doms_GtL_maps[i][j][3])(Eigen::all,Eigen::seq(Eigen::last,0,Eigen::fix<-1>));
        }
    }
    return Vi;
}


Matrix MultiDomain::global_to_gap_junct_left(const Matrix V) const
{
    Matrix Vg = Matrix::Zero(N_gap_junct_dofs,V.cols());
    for(auto i=0u;i<gap_junct_GtL_map.size();i++)
        Vg.block(gap_junct_GtL_map[i][2],0,gap_junct_GtL_map[i][3],V.cols())=
                V.block(gap_junct_GtL_map[i][1],0,gap_junct_GtL_map[i][3],V.cols());
    return Vg;
}

Matrix MultiDomain::global_to_gap_junct_right(const Matrix Vg) const
{
    Matrix V = Matrix::Zero(Vg.rows(),N_global_dofs);
    for(auto i=0u;i<gap_junct_GtL_map.size();i++)
        V.block(0,gap_junct_GtL_map[i][1],Vg.rows(),gap_junct_GtL_map[i][3])=
                Vg.block(0,gap_junct_GtL_map[i][2],Vg.rows(),gap_junct_GtL_map[i][3]);
    return V;
}

Matrix MultiDomain::gap_junct_to_global_left(const Matrix Vg) const
{
    Matrix V = Matrix::Zero(N_global_dofs,Vg.cols());
    for(auto i=0u;i<gap_junct_GtL_map.size();i++)
        V.block(gap_junct_GtL_map[i][1],0,gap_junct_GtL_map[i][3],Vg.cols())=
                Vg.block(gap_junct_GtL_map[i][2],0,gap_junct_GtL_map[i][3],Vg.cols());
    return V;
}

Matrix MultiDomain::gap_junct_to_global_right(const Matrix V) const
{
    Matrix Vg = Matrix::Zero(V.rows(),N_gap_junct_dofs);
    for(auto i=0u;i<gap_junct_GtL_map.size();i++)
        Vg.block(0,gap_junct_GtL_map[i][2],V.rows(),gap_junct_GtL_map[i][3])=
                V.block(0,gap_junct_GtL_map[i][1],V.rows(),gap_junct_GtL_map[i][3]);
    return Vg;
}

Vector MultiDomain::remove_measure_gap_junct(const Vector& Vg) const
{
    Vector Vg_nm = Vector::Zero(Vg.size());
    for(auto i=0u;i<gap_junct_GtL_map.size();i++)
        Vg_nm.segment(gap_junct_GtL_map[i][2],gap_junct_GtL_map[i][3])
         = Vg.segment(gap_junct_GtL_map[i][2],gap_junct_GtL_map[i][3])
                .cwiseQuotient(sides[gap_junct_GtL_map[i][0]].get_norm_dG()*sides[gap_junct_GtL_map[i][0]].get_dx());
    return Vg_nm;
}

Vector MultiDomain::add_measure_gap_junct(const Vector& Vg) const
{
    Vector Vg_nm = Vector::Zero(Vg.size());
    for(auto i=0u;i<gap_junct_GtL_map.size();i++)
        Vg_nm.segment(gap_junct_GtL_map[i][2],gap_junct_GtL_map[i][3])
         = Vg.segment(gap_junct_GtL_map[i][2],gap_junct_GtL_map[i][3])
                .cwiseProduct(sides[gap_junct_GtL_map[i][0]].get_norm_dG()*sides[gap_junct_GtL_map[i][0]].get_dx());
    return Vg_nm;
}

void MultiDomain::remove_free_constant(vector<Vector>& v_vec) const
{
    for(auto& v: v_vec)
        v = v.array()-v.array().mean();
}

vector<Matrix> MultiDomain::get_PS(vector<bool>& repeated) const
{
    vector<Matrix> PSp(N_domains);
    
    repeated.clear();
    repeated.resize(N_domains,false);

    if(N_domains>0)
        PSp[0] = doms[0].MyPoincareSteklow();

    for(auto i=1u;i<N_domains;i++)
    {
        if(doms[i-1]%doms[i])
        {
//            cout<<"Domains "<<i-1<<" and "<<i<<" are almost equal, up to a translation"<<flush<<endl;
            repeated[i] = true;
            PSp[i] = PSp[i-1];
        }
        else
            PSp[i] = doms[i].MyPoincareSteklow();
    }
    
    return PSp;
}

unsigned MultiDomain::get_N_domains() const
{
    return N_domains;
}

const Domain& MultiDomain::get_domain(unsigned i) const
{
    return doms[i];
}

Matrix MultiDomain::get_coords(unsigned i) const
{
    Matrix coords(N_local_dofs[i],2);
    
    for(auto j=0u;j<doms_GtL_maps[i].size();j++)
    {
        auto side = doms_GtL_maps[i][j];
        if(side[4]==1)
            coords.block(side[2],0,side[3],2) = sides[side[0]].get_G();
        else
            coords.block(side[2],0,side[3],2) = sides[side[0]].get_G()(Eigen::seq(Eigen::last,0,Eigen::fix<-1>),Eigen::all);
    }
    
    return coords;
}

Matrix MultiDomain::get_coords_gap_junct() const
{
    Matrix coords(N_gap_junct_dofs,2);
   
    for(auto j=0u;j<gap_junct_GtL_map.size();j++)
    {
        auto side = gap_junct_GtL_map[j];
        coords.block(side[2],0,side[3],2) = sides[side[0]].get_G();
    }
    
    return coords;
}

Vector MultiDomain::get_kappa_gap_junct() const
{
    Vector kappa(N_gap_junct_dofs);
   
    for(auto j=0u;j<gap_junct_GtL_map.size();j++)
    {
        auto side = gap_junct_GtL_map[j];
        kappa.segment(side[2],side[3]) = sides[side[0]].get_kappa();
    }
    
    return kappa;
}

Vector MultiDomain::get_angle_gap_junct() const
{
    Vector angle(N_gap_junct_dofs);
   
    for(auto j=0u;j<gap_junct_GtL_map.size();j++)
    {
        auto& side = gap_junct_GtL_map[j];
        unsigned n_side = sides[side[0]].get_n();
        const Matrix& N = sides[side[0]].get_N();
        Array angle_normal = N.block(0,0,n_side,1).array().acos();
        angle_normal = (N.block(0,1,n_side,1).array()>=0).select(angle_normal,-angle_normal);
        angle.segment(side[2],side[3]) 
                = (angle_normal>=0).select(angle_normal-M_PI/2.,angle_normal+M_PI/2.);
    }
    
    //we have angles from -pi/2 to pi/2
    
    return angle;
}


void MultiDomain::mesh()
{    
    for(auto& dom:doms)
        dom.mesh();
}

ostream& MultiDomain::write(ostream& os, string name)
{
    os<<name<<" = cell("<<N_domains<<",1);"<<endl;
    for(auto i=0u;i<N_domains;i++)
        doms[i].write(os,name+"{"+to_string(i+1)+"}");
    
    return os;
}

void MultiDomain::compute_boundary_errors(const vector<Vector>& traces, const vector<Vector>& derivs,
                                 const vector<Domain>& ref_doms, 
                                 const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                                 Array& err_traces, Array& err_derivs) const
{
    err_traces.resize(N_domains);
    err_derivs.resize(N_domains);
    for(unsigned i=0;i<N_domains;i++)
    {
        err_traces(i) = ref_doms[i].L2norm(ref_traces[i]-ref_doms[i].interpolate(traces[i],doms[i]));
        err_derivs(i) = ref_doms[i].L2norm(ref_derivs[i]-ref_doms[i].interpolate(derivs[i],doms[i]));
    }
}

void MultiDomain::compute_boundary_errors_inner(const vector<Vector>& traces, const vector<Vector>& derivs,
                                 const vector<Domain>& ref_doms, 
                                 const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                                 Array& err_traces, Array& err_derivs) const
{
    err_traces.resize(N_domains);
    err_derivs.resize(N_domains);
    for(unsigned i=0;i<N_domains;i++)
    {
        err_traces(i) = ref_doms[i].L2norm_inner(ref_traces[i]-ref_doms[i].interpolate_inner(traces[i],doms[i]));
        err_derivs(i) = ref_doms[i].L2norm_inner(ref_derivs[i]-ref_doms[i].interpolate_inner(derivs[i],doms[i]));
    }
}

void MultiDomain::compute_boundary_errors(const vector<Vector>& traces, const vector<Vector>& derivs,
                                 const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                                 Array& err_traces, Array& err_derivs) const
{
    err_traces.resize(N_domains);
    err_derivs.resize(N_domains);
    for(unsigned i=0;i<N_domains;i++)
    {
        err_traces(i) = doms[i].L2norm(ref_traces[i]-traces[i]);
        err_derivs(i) = doms[i].L2norm(ref_derivs[i]-derivs[i]);
    }
}

void MultiDomain::compute_boundary_errors_inner(const vector<Vector>& traces, const vector<Vector>& derivs,
                                 const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                                 Array& err_traces, Array& err_derivs) const
{
    err_traces.resize(N_domains);
    err_derivs.resize(N_domains);
    for(unsigned i=0;i<N_domains;i++)
    {
        err_traces(i) = doms[i].L2norm_inner(ref_traces[i]-traces[i]);
        err_derivs(i) = doms[i].L2norm_inner(ref_derivs[i]-derivs[i]);
    }
}

Real MultiDomain::compute_domain_error(const Domain& dom, const Vector& trace, const Vector& deriv,
                                       const Domain& ref_dom, const Vector& ref_trace, const Vector& ref_deriv) const
{
    vector<Real> points = dom.get_non_boundary_points();

    vector<Real> u = dom.eval_on_points(points,trace,deriv);
    vector<Real> ref_u = ref_dom.eval_on_points(points,ref_trace,ref_deriv);
    
    Real err=0.;
    Real ref_norm=0.;
    for(unsigned i=0;i<u.size();i++)
    {
        err += pow(u[i]-ref_u[i],2);
        ref_norm += pow(ref_u[i],2);
    }
    ref_norm = sqrt(ref_norm);
    err = sqrt(err);
    
    Real tol = 1e-10;
    if(ref_norm>tol)
        return err/ref_norm;
    else
        return err;
}

void MultiDomain::compute_domain_errors(const vector<Vector>& traces, const vector<Vector>& derivs,
                               const vector<Domain>& ref_doms, 
                               const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                               Array& err_doms) const
{    
    err_doms.resize(N_domains);
    for(unsigned i=0;i<N_domains;i++)
        err_doms[i] = compute_domain_error(doms[i],traces[i],derivs[i],
                                           ref_doms[i],ref_traces[i],ref_derivs[i]);          
}

void MultiDomain::compute_domain_errors(const vector<Vector>& traces, const vector<Vector>& derivs,
                               const vector<Vector>& ref_traces, const vector<Vector>& ref_derivs,
                               Array& err_doms) const
{    
    err_doms.resize(N_domains);
    for(unsigned i=0;i<N_domains;i++)
        err_doms[i] = compute_domain_error(doms[i],traces[i],derivs[i],
                                           doms[i],ref_traces[i],ref_derivs[i]);          
}

unsigned MultiDomain::get_n_dofs() const
{
    return N_global_dofs;
}

unsigned MultiDomain::get_n_transmembrane_dofs() const
{
    return N_transmembrane_dofs;
}

unsigned MultiDomain::get_n_outer_bnd_dofs() const
{
    return N_outer_bnd_dofs;
}

void MultiDomain::save_domains(const string ref_name) const
{
    ofstream ofile(ref_name+"_domains.bin", std::ios::out | std::ios::binary | std::ios::trunc);    
    for(auto& dom:doms)
        dom.save(ofile);    
    ofile.close();
}

void MultiDomain::save_boundary_data(const string ref_name, const vector<Vector>& traces, const vector<Vector>& derivs) const
{
    ofstream ofile(ref_name+"_solutions.bin", std::ios::out | std::ios::binary | std::ios::trunc);  
    for(unsigned i=0;i<N_domains;i++)
    {
        write_binary(ofile, traces[i]);
        write_binary(ofile, derivs[i]);
    }
    ofile.close();
}

void MultiDomain::save(const string ref_name, const vector<Vector>& traces, const vector<Vector>& derivs) const
{
    save_domains(ref_name);
    save_boundary_data(ref_name,traces,derivs);
}

void MultiDomain::load(const string ref_name, vector<Vector>& traces, vector<Vector>& derivs,  vector<Domain>& domains) const
{
    domains = load_domains(ref_name);
    load_boundary_data(ref_name,traces,derivs);
}

void MultiDomain::get_multidomain_desc(const int n, Parameters& param,
                  Boundary& outer_bnd, vector<Side>& sides, vector<SidesToBoundaryMap>& inner_dom_desc)
{
    /*  List of geometries. 
     * For n<=0 cell size is not in same as scale as real cells.
     * For n>=1 we consider the same scale and more involved problems.
     * 
     * n = -7: One cell made of two half circles connected by two segments. Outer boundary is an ellipse.
     * n = -6: One ellipse as inner cell and another as outer boundary.
     * n = -5: Two non touching inner domains, smoothed half circles. Outer domain is a circle.
     * n = -4: Two non touching inner domains, half circles. Outer domain is a circle.
     * n = -3: Two touching inner domains, half circles. Outer domain is a circle.
     * n = -2: One inner domain, a circles. Outer domain is a circle.
     * n = -1: Four touching inner cells. Outer domain is an ellipse.
     * n =  0: Two clusters of two cells each. Outer domain is an ellipse.
     * n =  1: Example given by Fatemeh but with different scale. It is a spiral of cells. 
     *         For our purpose here we use a smaller scale, i.e. the same as real cells, to see stiffness.
     * n =  2: One cluster of ten cells having shape similar to real ones.
     */
    
    
//    if(n==-7)
//    {
//        // Outer Domain = ellipse 
//        // Inner domain = two half circles connected by segments
//        Real x = 1.;
//        Real r = 1.5;
//        Boundary outer_bnd = Boundary::get_ellipse({0.,0.},4.,3.,dG,1.);
//        Side s0 = Arc({-x,0.},r,M_PI/2.,1.5*M_PI,dG,1.);
//        Side s1 = Segment({-x,-r},{x,-r},dG,1.);
//        Side s2 = Arc({x,0.},r,1.5*M_PI,2.5*M_PI,dG,1.);
//        Side s3 = Segment({x,r},{-x,r},dG,1.);
//        SidesToBoundaryMap inner_dom_1(4);
//        inner_dom_1[0]={0,1};
//        inner_dom_1[1]={1,1};
//        inner_dom_1[2]={2,1};
//        inner_dom_1[3]={3,1};
//        vector<Side> sides={s0,s1,s2,s3};
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1});
//        return mdom;
//    }
//    
//    if(n==-6)
//    {
//        Boundary outer_bnd = Boundary::get_ellipse({0.,0.},4.,3.,dG,1.);
//        Side s0 = EllipseArc({0.,0.},2.,0.5,0.,2*M_PI,dG,1.);
//        SidesToBoundaryMap inner_dom_1(1);
//        inner_dom_1[0]={0,1};
//        vector<Side> sides={s0};
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1});
//        return mdom;
//    }
//    
//    if(n==-5)
//    {
//        // Outer domain = circle radius 4
//        // Two non touching inner domains: smoothed half circles of radius 2
//        Boundary outer_bnd = Boundary::get_circle({0.,0.},4.,dG,1.);
//        Point p{0.2,0.};
//        Real r = 0.2;
//        Side s0 = -p+Arc({-r,0.},2.,0.5*M_PI,1.5*M_PI,dG,1.);
//        Side s1 = -p+Arc({-r,-2.+r},r,-0.5*M_PI,0.,dG,1.);
//        Side s2 = -p+Segment({0.,-2.+r},{0.,2.-r},dG,1.);
//        Side s3 = -p+Arc({-r,2.-r},r,0.,0.5*M_PI,dG,1.);
//        
//        Side s4 = p+Arc({r,0.},2.,-0.5*M_PI,0.5*M_PI,dG,1.);
//        Side s5 = p+Arc({r,2.-r},r,0.5*M_PI,M_PI,dG,1.);
//        Side s6 = p+Segment({0.,2.-r},{0.,-2.+r},dG,1.);       
//        Side s7 = p+Arc({r,-2.+r},r,M_PI,1.5*M_PI,dG,1.);
//        SidesToBoundaryMap inner_dom_1(4);
//        inner_dom_1[0]={0,1};
//        inner_dom_1[1]={1,1};
//        inner_dom_1[2]={2,1};
//        inner_dom_1[3]={3,1};
//        SidesToBoundaryMap inner_dom_2(4);
//        inner_dom_2[0]={4,1};
//        inner_dom_2[1]={5,1};
//        inner_dom_2[2]={6,1};
//        inner_dom_2[3]={7,1};
//        vector<Side> sides={s0,s1,s2,s3,s4,s5,s6,s7};
//        
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1,inner_dom_2});
//        return mdom;
//    }
//    
//    if(n==-4)
//    {
//        // Outer domain = circle radius 4
//        // Two non touching inner domains: half circle of radius 2
//        Boundary outer_bnd = Boundary::get_circle({0.,0.},4.,dG,1.);
//        Point p{0.2,0.};
//        Side s0 = -p+Arc({0.,0.},2.,0.5*M_PI,1.5*M_PI,dG,1.);
//        Side s1 = -p+Segment({0.,-2.},{0.,2.},dG,1.);
//        Side s2 = p+Arc({0.,0.},2.,-0.5*M_PI,0.5*M_PI,dG,1.);
//        Side s3 = p+Segment({0.,2.},{0.,-2.},dG,1.);       
//        SidesToBoundaryMap inner_dom_1(2);
//        inner_dom_1[0]={0,1};
//        inner_dom_1[1]={1,1};
//        SidesToBoundaryMap inner_dom_2(2);
//        inner_dom_2[0]={2,1};
//        inner_dom_2[1]={3,1};
//        vector<Side> sides={s0,s1,s2,s3};
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1,inner_dom_2});
//        return mdom;
//    }
//    
//    if(n==-3)
//    {
//        // Outer domain = circle radius 4
//        // Two touching inner domains: half circles of radius 2
//        Boundary outer_bnd = Boundary::get_circle({0.,0.},4.,dG,1.);
//        Side s0 = Arc({0.,0.},2.,0.5*M_PI,1.5*M_PI,dG,1.);       
//        Side s1 = Segment({0.,-2.},{0.,2.},dG,1.);
//        Side s2 = Arc({0.,0.},2.,-0.5*M_PI,0.5*M_PI,dG,1.);
//        SidesToBoundaryMap inner_dom_1(2);
//        inner_dom_1[0]={0,1};
//        inner_dom_1[1]={1,1};
//        SidesToBoundaryMap inner_dom_2(2);
//        inner_dom_2[0]={1,0};
//        inner_dom_2[1]={2,1};
//        vector<Side> sides={s0,s1,s2};
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1,inner_dom_2});
//        return mdom;
//    }
//    
//    if(n==-2)
//    {
//        // Outer Domain = circle radius 4
//        // Inner domain = circle radius 2
//        Boundary outer_bnd = Boundary::get_circle({0.,0.},4.,dG,1.);
//        Side s0 = Arc({0.,0.},2.,0.,2*M_PI,dG,1.);
//        SidesToBoundaryMap inner_dom_1(1);
//        inner_dom_1[0]={0,1};
//        vector<Side> sides={s0};
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1});
//        return mdom;
//    }
//    
//    if(n==-1)
//    {
//        // Outer domain = ellipse
//        // four touching cells
//        Real r = 1.2;                
//        Boundary outer_bnd = Boundary::get_ellipse({0.,0.},4.*r,2.*r,dG,1.);
//        
//        Side q1 = Segment({0.,-r},{0.,r},dG,1.);//0
//        Side q2 = Segment({0.,r},{-2.*r,r},dG,1.);//1
//        Side q3 = Segment({-2.*r,r},{-2.*r,-r},dG,1.);//2
//        Side q4 = Segment({-2.*r,-r},{0.,-r},dG,1.);//3
//        Side a0 = Arc({-2*r,0.},r,0.5*M_PI,1.5*M_PI,dG,1.);//4        
//        
//        Side p1 = Segment({0.,-r},{2.*r,-r},dG,1.);//5
//        Side p2 = Segment({2.*r,-r},{2.*r,r},dG,1.);//6
//        Side p3 = Segment({2.*r,r},{0.,r},dG,1.);//7
//        //0 reversed
//        Side a1 = Arc({2.*r,0.},r,-0.5*M_PI,0.5*M_PI,dG,1.);//8
//        
//        SidesToBoundaryMap inner_dom_1(4);
//        inner_dom_1[0]={0,1};
//        inner_dom_1[1]={1,1};
//        inner_dom_1[2]={2,1};
//        inner_dom_1[3]={3,1};
//        SidesToBoundaryMap inner_dom_2(4);
//        inner_dom_2[0]={5,1};
//        inner_dom_2[1]={6,1};
//        inner_dom_2[2]={7,1};
//        inner_dom_2[3]={0,0};
//        SidesToBoundaryMap inner_dom_3(2);
//        inner_dom_3[0]={2,0};
//        inner_dom_3[1]={4,1};
//        SidesToBoundaryMap inner_dom_4(2);
//        inner_dom_4[0]={6,0};
//        inner_dom_4[1]={8,1};
//        vector<Side> sides={q1,q2,q3,q4,a0,p1,p2,p3,a1};
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1,inner_dom_2,inner_dom_3,inner_dom_4});
//        return mdom;
//    }
//    
//    if(n==0)
//    {
//        //same cells as in -1 but now separeted in two clusters
//        // Outer domain = ellipse
//        // two clusters of two cells
//        Real r = 1.2;
//        Real d1 = 0.68;
//        Real d2 = 0.*4./3.;
//        Point p{d1/2.,d2/2.};
//        
//        Boundary outer_bnd = Boundary::get_ellipse({0.,0.},4.*r+d1,2.*r+d2,dG,1.);
//        
//        Side q1 = -p+Segment({0.,-r},{0.,r},dG,1.);//0
//        Side q2 = -p+Segment({0.,r},{-2.*r,r},dG,1.);//1
//        Side q3 = -p+Segment({-2.*r,r},{-2.*r,-r},dG,1.);//2
//        Side q4 = -p+Segment({-2.*r,-r},{0.,-r},dG,1.);//3
//        Side a0 = -p+Arc({-2*r,0.},r,0.5*M_PI,1.5*M_PI,dG,1.);//4        
//        
//        Side p1 = p+Segment({0.,-r},{2.*r,-r},dG,1.);//5
//        Side p2 = p+Segment({2.*r,-r},{2.*r,r},dG,1.);//6
//        Side p3 = p+Segment({2.*r,r},{0.,r},dG,1.);//7
//        Side p4 = p+Segment({0.,r},{0.,-r},dG,1.);//8
//        Side a1 = p+Arc({2.*r,0.},r,-0.5*M_PI,0.5*M_PI,dG,1.);//9
//        
//        SidesToBoundaryMap inner_dom_1(4);
//        inner_dom_1[0]={0,1};
//        inner_dom_1[1]={1,1};
//        inner_dom_1[2]={2,1};
//        inner_dom_1[3]={3,1};
//        SidesToBoundaryMap inner_dom_2(4);
//        inner_dom_2[0]={5,1};
//        inner_dom_2[1]={6,1};
//        inner_dom_2[2]={7,1};
//        inner_dom_2[3]={8,1};
//        SidesToBoundaryMap inner_dom_3(2);
//        inner_dom_3[0]={2,0};
//        inner_dom_3[1]={4,1};
//        SidesToBoundaryMap inner_dom_4(2);
//        inner_dom_4[0]={6,0};
//        inner_dom_4[1]={9,1};
//        vector<Side> sides={q1,q2,q3,q4,a0,p1,p2,p3,p4,a1};
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1,inner_dom_2,inner_dom_3,inner_dom_4});
//        return mdom;
//    }
//    
//    
//    if(n==1)
//    {
//        // Use u = 0.125 to recover Fatemeh example
//        // We chose scale so that cells have length 
//        // ~100 um = 100e-6 m = 10000e-6 cm = 1e-2 cm
////        Real u = 0.125; //cm.
//        Real u = 5e-3; //cm.
//        Boundary outer_bnd = Boundary::get_square({4.*u,4.*u},8.*u,dG,1.);
//        
//        // cells have side of length 2u=1e-2cm
//        Segment h1 = Segment({0,0},{2.*u,0.},dG,1.);
//        Segment h2 = Segment({0,0},{-2.*u,0.},dG,1.);
//        Segment v1 = Segment({0,0},{0.,2.*u},dG,1.);
//        Segment v2 = Segment({0,0},{0.,-2.*u},dG,1.);
//        Segment d1 = Segment({0,0},{u,u},dG,1.);
//        Segment d2 = Segment({0,0},{-u,-u},dG,1.);
//                
//        Side s0 = u*Point{1.,1.}+d1;
//        Side s1 = u*Point{2.,2.}+v1;
//        Side s2 = u*Point{2.,4.}+d2;
//        Side s3 = u*Point{1.,3.}+v2;
//        
//        SidesToBoundaryMap inner_dom_1(4);
//        inner_dom_1[0]={0,1};
//        inner_dom_1[1]={1,1};
//        inner_dom_1[2]={2,1};
//        inner_dom_1[3]={3,1};
//        
//        Side s4 = u*Point{1.,1.}+h1;
//        Side s5 = u*Point{3.,1.}+d1;
//        Side s6 = u*Point{4.,2.}+h2;
//        
//        SidesToBoundaryMap inner_dom_2(4);
//        inner_dom_2[0]={4,1};
//        inner_dom_2[1]={5,1};
//        inner_dom_2[2]={6,1};
//        inner_dom_2[3]={0,0};
//        
//        Side s7 = u*Point{3.,1.}+h1;
//        Side s8 = u*Point{5.,1.}+d1;
//        Side s9 = u*Point{6.,2.}+h2;
//        
//        SidesToBoundaryMap inner_dom_3(4);
//        inner_dom_3[0]={7,1};
//        inner_dom_3[1]={8,1};
//        inner_dom_3[2]={9,1};
//        inner_dom_3[3]={5,0};
//        
//        Side s10 = u*Point{5.,1.}+h1;
//        Side s11 = u*Point{7.,1.}+v1;
//        Side s12 = u*Point{7.,3.}+d2;
//        
//        SidesToBoundaryMap inner_dom_4(4);
//        inner_dom_4[0]={10,1};
//        inner_dom_4[1]={11,1};
//        inner_dom_4[2]={12,1};
//        inner_dom_4[3]={8,0};
//        
//        Side s13 = u*Point{7.,3.}+v1;
//        Side s14 = u*Point{7.,5.}+d2;
//        Side s15 = u*Point{6.,4.}+v2;
//        
//        SidesToBoundaryMap inner_dom_5(4);
//        inner_dom_5[0]={12,0};
//        inner_dom_5[1]={13,1};
//        inner_dom_5[2]={14,1};
//        inner_dom_5[3]={15,1};
//        
//        Side s16 = u*Point{7.,5.}+v1;
//        Side s17 = u*Point{7.,7.}+d2;
//        Side s18 = u*Point{6.,6.}+v2;
//        
//        SidesToBoundaryMap inner_dom_6(4);
//        inner_dom_6[0]={14,0};
//        inner_dom_6[1]={16,1};
//        inner_dom_6[2]={17,1};
//        inner_dom_6[3]={18,1};
//        
//        Side s19 = u*Point{7.,7.}+h2;
//        Side s20 = u*Point{5.,7.}+d2;
//        Side s21 = u*Point{4.,6.}+h1;
//        
//        SidesToBoundaryMap inner_dom_7(4);
//        inner_dom_7[0]={17,0};
//        inner_dom_7[1]={19,1};
//        inner_dom_7[2]={20,1};
//        inner_dom_7[3]={21,1};
//        
//        Side s22 = u*Point{5.,7.}+h2;
//        Side s23 = u*Point{3.,7.}+d2;
//        Side s24 = u*Point{2.,6.}+h1;
//        
//        SidesToBoundaryMap inner_dom_8(4);
//        inner_dom_8[0]={20,0};
//        inner_dom_8[1]={22,1};
//        inner_dom_8[2]={23,1};
//        inner_dom_8[3]={24,1};
//        
//        Side s25 = u*Point{3.,7.}+h2;
//        Side s26 = u*Point{1.,7.}+v2;
//        Side s27 = u*Point{1.,5.}+d1;
//        
//        SidesToBoundaryMap inner_dom_9(4);
//        inner_dom_9[0]={23,0};
//        inner_dom_9[1]={25,1};
//        inner_dom_9[2]={26,1};
//        inner_dom_9[3]={27,1};
//        
//        vector<Side> sides={s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27};
//        MultiDomain mdom(outer_bnd,sides,{inner_dom_1,inner_dom_2,inner_dom_3,inner_dom_4,inner_dom_5,inner_dom_6,inner_dom_7,inner_dom_8,inner_dom_9});
//        return mdom;
//        
//    }
//    
//    if(n==2)
//    {
//        // One cluster of ten cells with complicated shape.
//        Real u = 0.125;
////        Real u = 5e-3;
////        Boundary outer_bnd = Boundary::get_ellipse({13*u,1*u},32.*u,16.*u,dG);
//        Boundary outer_bnd = Boundary::get_rectangle({-7*u,-8.6*u},{33*u,10.6*u},dG,1.);
//        
//        auto seg = [=](Real x1, Real y1, Real x2, Real y2){return Segment({x1*u,y1*u},{x2*u,y2*u},dG,1.);};
//        vector<Side> sides;
//        
//        Real r = 0.25;
//        Point c1{5+r*(2.+sqrt(5.)),1+r};
//        Point c2{19+r*(2.+sqrt(5.)),4+r};
//        Real alpha = M_PI/2.+acos(2./sqrt(5.));
//        Real beta = 1.5*M_PI;
//        
//        Point c3{18+r*(1.+sqrt(2.)),-1-r};
//        Real gamma = M_PI/2.;
//        Real delta = M_PI+M_PI/4.;//acos(-1./sqrt(2.));
//                
//        sides.push_back(seg(-3,0,0,0)); // this segment is 150um if u=5e-3
//        sides.push_back(seg(0,0,1,2));
//        sides.push_back(seg(1,2,-2,2));
//        sides.push_back(seg(-2,2,-3,0));
//        sides.push_back(seg(0,0,6,-1));
//        sides.push_back(seg(6,-1,7,1));
//        sides.push_back(seg(7,1,c1(0),1));//6
//        sides.push_back(-Arc(c1*u,r*u,alpha,beta,dG,1.));
//        sides.push_back(seg(c1(0)-r/sqrt(5),c1(1)+2.*r/sqrt(5.),9,3));
//        sides.push_back(seg(9,3,10,5));//9
//        sides.push_back(seg(10,5,7,5));
//        sides.push_back(seg(7,5,1,2));
//        sides.push_back(seg(9,3,15,3));
//        sides.push_back(seg(15,3,16,5));
//        sides.push_back(seg(16,5,10,5));
//        sides.push_back(seg(15,3,17,2));//15
//        sides.push_back(seg(17,2,20,2)); 
//        sides.push_back(seg(20,2,21,4));
//        sides.push_back(seg(21,4,c2(0),4));//18
//        sides.push_back(-Arc(c2*u,r*u,alpha,beta,dG,1.));
//        sides.push_back(seg(c2(0)-r/sqrt(5),c2(1)+2.*r/sqrt(5.),23,6));
//        sides.push_back(seg(23,6,24,8));
//        sides.push_back(seg(24,8,22,8));
//        sides.push_back(seg(22,8,16,5));//23
//        sides.push_back(seg(23,6,27,6));
//        sides.push_back(seg(27,6,28,8));
//        sides.push_back(seg(28,8,24,8));//26
//        sides.push_back(seg(20,2,26,2));
//        sides.push_back(seg(26,2,27,4));
//        sides.push_back(seg(27,4,21,4));//29
//        sides.push_back(seg(6,-1,12,-1));
//        sides.push_back(seg(12,-1,13,1));
//        sides.push_back(seg(13,1,7,1));//32
//        sides.push_back(seg(12,-1,16,-1));
//        sides.push_back(seg(16,-1,18,-3));
//        sides.push_back(seg(18,-3,20,-3));
//        sides.push_back(seg(20,-3,c3(0)-r/sqrt(2.),c3(1)-r/sqrt(2.)));
//        sides.push_back(-Arc(c3*u,r*u,gamma,delta,dG,1.));
//        sides.push_back(seg(c3(0),c3(1)+r,20,-1));
//        sides.push_back(seg(20,-1,21,1));
//        sides.push_back(seg(21,1,13,1));//40
//        sides.push_back(seg(18,-3,15,-3));
//        sides.push_back(seg(15,-3,14,-5));
//        sides.push_back(seg(14,-5,21,-5));
//        sides.push_back(seg(21,-5,22,-3));
//        sides.push_back(seg(22,-3,20,-3));//45
//        sides.push_back(seg(21,-5,25,-5));
//        sides.push_back(seg(25,-5,26,-3));
//        sides.push_back(seg(26,-3,22,-3));//48
//        
//        
//        vector<SidesToBoundaryMap> inner_doms;
//        inner_doms.push_back(SidesToBoundaryMap(4));
//        inner_doms.back()[0]={0,1};
//        inner_doms.back()[1]={1,1};
//        inner_doms.back()[2]={2,1};
//        inner_doms.back()[3]={3,1};
//        inner_doms.push_back(SidesToBoundaryMap(9));
//        inner_doms.back()[0]={1,0};
//        inner_doms.back()[1]={4,1};
//        inner_doms.back()[2]={5,1};
//        inner_doms.back()[3]={6,1};
//        inner_doms.back()[4]={7,1};
//        inner_doms.back()[5]={8,1};
//        inner_doms.back()[6]={9,1};
//        inner_doms.back()[7]={10,1};
//        inner_doms.back()[8]={11,1};
//        inner_doms.push_back(SidesToBoundaryMap(4));
//        inner_doms.back()[0]={9,0};
//        inner_doms.back()[1]={12,1};
//        inner_doms.back()[2]={13,1};
//        inner_doms.back()[3]={14,1};
//        inner_doms.push_back(SidesToBoundaryMap(10));
//        inner_doms.back()[0]={13,0};
//        inner_doms.back()[1]={15,1};
//        inner_doms.back()[2]={16,1};
//        inner_doms.back()[3]={17,1};
//        inner_doms.back()[4]={18,1};
//        inner_doms.back()[5]={19,1};
//        inner_doms.back()[6]={20,1};
//        inner_doms.back()[7]={21,1};
//        inner_doms.back()[8]={22,1};
//        inner_doms.back()[9]={23,1};
//        inner_doms.push_back(SidesToBoundaryMap(4));
//        inner_doms.back()[0]={21,0};
//        inner_doms.back()[1]={24,1};
//        inner_doms.back()[2]={25,1};
//        inner_doms.back()[3]={26,1};
//        inner_doms.push_back(SidesToBoundaryMap(4));
//        inner_doms.back()[0]={17,0};
//        inner_doms.back()[1]={27,1};
//        inner_doms.back()[2]={28,1};
//        inner_doms.back()[3]={29,1};
//        inner_doms.push_back(SidesToBoundaryMap(4));
//        inner_doms.back()[0]={5,0};
//        inner_doms.back()[1]={30,1};
//        inner_doms.back()[2]={31,1};
//        inner_doms.back()[3]={32,1};
//        inner_doms.push_back(SidesToBoundaryMap(9));
//        inner_doms.back()[0]={31,0};
//        inner_doms.back()[1]={33,1};
//        inner_doms.back()[2]={34,1};
//        inner_doms.back()[3]={35,1};
//        inner_doms.back()[4]={36,1};
//        inner_doms.back()[5]={37,1};
//        inner_doms.back()[6]={38,1};
//        inner_doms.back()[7]={39,1};
//        inner_doms.back()[8]={40,1};
//        inner_doms.push_back(SidesToBoundaryMap(6));
//        inner_doms.back()[0]={35,0};
//        inner_doms.back()[1]={41,1};
//        inner_doms.back()[2]={42,1};
//        inner_doms.back()[3]={43,1};
//        inner_doms.back()[4]={44,1};
//        inner_doms.back()[5]={45,1};
//        inner_doms.push_back(SidesToBoundaryMap(4));
//        inner_doms.back()[0]={44,0};
//        inner_doms.back()[1]={46,1};
//        inner_doms.back()[2]={47,1};
//        inner_doms.back()[3]={48,1};
//        
//        MultiDomain mdom(outer_bnd,sides,inner_doms);
//        return mdom;        
//    }
    
    if(n==3)
    {
        Real cl = 1.0*1e-2;
        Real cw = 0.20*1e-2;        
        Point p0(0.,0.); //origin of first cell of first cluster
        Point d1(1e-4,0.); //distance 1 between clusters
        Point d2(0.,cw/2.); //distance 2 between clusters
        unsigned nx_clusters = 1;
        unsigned ny_clusters = 1;
        unsigned nx_cells_per_cluster = 1;
        unsigned ny_cells_per_cluster = 20;   
        Real dG = param.P20_dG;
        Real dGl = dG;
        Real dGw = dG;
        Real vertamp = 0.5e-4;
        unsigned vertfreq = 0.;        
        Real horamp = 0.;
        Real horlen = 20e-4;
        Real horpos = 0.5;
        int horalt = -1;
            
        if(param.P20_cell_l>0.)
            cl=param.P20_cell_l;
        else
            param.P20_cell_l=cl;
        if(param.P20_cell_w>0.)
            cw=param.P20_cell_w;
        else
            param.P20_cell_w=cw;
        if(param.P20_nx>0.)
            nx_cells_per_cluster=param.P20_nx;
        else
            param.P20_nx=nx_cells_per_cluster;
        if(param.P20_ny>0.)
            ny_cells_per_cluster=param.P20_ny;
        else
            param.P20_ny=ny_cells_per_cluster;
        if(param.P20_vertamp>0.)
            vertamp=param.P20_vertamp;
        else
            param.P20_vertamp=vertamp;
        if(param.P20_vertfreq>0.)
            vertfreq=param.P20_vertfreq;
        else
            param.P20_vertfreq=vertfreq;
        if(param.P20_horamp>0.)
            horamp=param.P20_horamp;
        else
            param.P20_horamp=horamp;
        if(param.P20_horlen>0.)
            horlen=param.P20_horlen;
        else
            param.P20_horlen=horlen;
        if(param.P20_horpos>0.)
            horpos=param.P20_horpos;
        else
            param.P20_horpos=horpos;
        if(param.P20_horalt>-2)
            horalt=param.P20_horalt;
        else
            param.P20_horalt=horalt;
        
        cout<<"Building Sub Domains..."<<endl;
        
        vector<SidesToBoundaryMap> cells(0);
        sides.resize(0);
        
        Point p1;
        for(unsigned j=0;j<ny_clusters;j++)
        {
            p1 = p0 +Point(0.,j*cw*ny_cells_per_cluster) +j*d2 ;
            for(unsigned i=0;i<nx_clusters;i++)
            { 
                add_cells_cluster(param,p1,cells,sides);
//                add_cells_cluster_with_holes(cl,cw,pl,pw,p1,nx_cells_per_cluster,ny_cells_per_cluster,
//                                  dGl,dGw,cells,sides);
                p1 += Point(nx_cells_per_cluster*cl,0.)+d1;
            }
        }
        p1 = p0 +Point(nx_clusters*cl*nx_cells_per_cluster,ny_clusters*cw*ny_cells_per_cluster)
                +(nx_clusters-1)*d1+(ny_clusters-1)*d2;
                                        
        // Definition of the outer boundary
        Real b1 = 10.*cl;
        Real b2 = 10.*cw;
        Point q0 = p0 - Point(b1,b2);
        Point q1 = p1 + Point(b1,b2);
        Real outer_dG = 10.*dG;
        Segment s1(q0,{q1(0),q0(1)},outer_dG,1.,0.);
        Segment s2({q1(0),q0(1)},q1,outer_dG,1.,0.);
        Segment s3(q1,{q0(0),q1(1)},outer_dG,1.,0.);
        Segment s4({q0(0),q1(1)},q0,outer_dG,1.,0.);
        outer_bnd = Boundary({s1,s2,s3,s4});
        
        inner_dom_desc=cells;
        return;

    }
    if(n==4)
    {
        Real cl = 1e-4;
        Real h = 1./3.; 
        Real dG = param.P20_dG;       
        Real kappa = 1./param.P20_Rl;
            
        cout<<"Building Sub Domains..."<<endl;
        
        vector<SidesToBoundaryMap> cells(0);
        sides.resize(0);
        
        cells.resize(0);
        sides.resize(0);
        
        auto seg = [=](Real x1, Real y1, Real x2, Real y2)
                    {return Segment({x1*cl,(723.-y1)*cl},{x2*cl,(723.-y2)*cl},dG,1.,kappa);};
        
        sides.push_back( seg(3,720,93,720) );//0
        sides.push_back( seg(93,720,27,636) );//1
        sides.push_back( seg(27,636,3,720) );//2
        cells.push_back(SidesToBoundaryMap(3));
        cells.back()[0] = {0,1};
        cells.back()[1] = {1,1};
        cells.back()[2] = {2,1};
        
        sides.push_back( seg(93,720,153,440) );//3
        sides.push_back( seg(153,440,90,347) );//4
        sides.push_back( seg(90,347,27,636) );//5
        cells.push_back(SidesToBoundaryMap(4));
        cells.back()[0] = {1,0};
        cells.back()[1] = {3,1};
        cells.back()[2] = {4,1};
        cells.back()[3] = {5,1};
        
        sides.push_back( seg(153,440,177,336) );//6
        sides.push_back( seg(177,336,238,0) );//7
        sides.push_back( seg(238,0,156,0) );//8
        sides.push_back( seg(156,0,90,347) );//9
        cells.push_back(SidesToBoundaryMap(5));
        cells.back()[0] = {4,0};
        cells.back()[1] = {6,1};
        cells.back()[2] = {7,1};
        cells.back()[3] = {8,1};
        cells.back()[4] = {9,1};
        
        sides.push_back( seg(153,440,210,490) );//10
        sides.push_back( seg(210,490,280,505) );//11
        sides.push_back( seg(280,505,354,481) );//12
        sides.push_back( seg(354,481,411,441) );//13
        sides.push_back( seg(411,441,403,378) );//14 --
        sides.push_back( seg(403,378,343,410) );//15
        sides.push_back( seg(343,410,265,427) );//16
        sides.push_back( seg(265,427,177,336) );//17
        cells.push_back(SidesToBoundaryMap(9));
        cells.back()[0] = {6,0};
        cells.back()[1] = {10,1};
        cells.back()[2] = {11,1};
        cells.back()[3] = {12,1};
        cells.back()[4] = {13,1};
        cells.back()[5] = {14,1};
        cells.back()[6] = {15,1};
        cells.back()[7] = {16,1};
        cells.back()[8] = {17,1};        
        
        
        sides.push_back( seg(411,441,480,505) );//18
        sides.push_back( seg(480,505,552,494) );//19
        sides.push_back( seg(552,494,565,417) );//20
        sides.push_back( seg(565,417,497,432) );//21
        sides.push_back( seg(497,432,498,357) );//22
        sides.push_back( seg(498,357,523,217) );//23
        sides.push_back( seg(523,217,535,165) );//24
        sides.push_back( seg(535,165,552,82) );//25 --
        sides.push_back( seg(552,82,570,0) );//26
        sides.push_back( seg(570,0,487,0) );//27 sopra
        sides.push_back( seg(487,0,457,165) );//28
        sides.push_back( seg(457,165,434,267) );//29
        sides.push_back( seg(434,267,403,378) );//30
        cells.push_back(SidesToBoundaryMap(14));
        cells.back()[0] = {14,0};
        cells.back()[1] = {18,1};
        cells.back()[2] = {19,1};
        cells.back()[3] = {20,1};
        cells.back()[4] = {21,1};
        cells.back()[5] = {22,1};
        cells.back()[6] = {23,1};
        cells.back()[7] = {24,1};
        cells.back()[8] = {25,1};
        cells.back()[9] = {26,1};
        cells.back()[10] = {27,1};
        cells.back()[11] = {28,1};
        cells.back()[12] = {29,1};
        cells.back()[13] = {30,1};
        
        //connessione 
        Real dx = -80;
        Real dy = 0;
        sides.push_back( seg(535,165,639+dx,165+dy) );//31
        sides.push_back( seg(639+dx,165+dy,658+dx,82+dy) );//32-
        sides.push_back( seg(658+dx,82+dy,552,82) );//33        
        cells.push_back(SidesToBoundaryMap(4));
        cells.back()[0] = {25,0};
        cells.back()[1] = {31,1};
        cells.back()[2] = {32,1};
        cells.back()[3] = {33,1};
//        
//        // cuore
        Point p0(dx*cl,-dy*cl);
        sides.push_back( p0 + seg(639,165,649,249) );//34
        sides.push_back( p0 + seg(649,249,730,372) );//35
        sides.push_back( p0 + seg(730,372,825,456) );//36
        sides.push_back( p0 + seg(825,456,915,514) );//37
        sides.push_back( p0 + seg(915,514,1057,420) );//38
        sides.push_back( p0 + seg(1057,420,1134,343) );//39
        sides.push_back( p0 + seg(1134,343,1195,234) );//40
        sides.push_back( p0 + seg(1195,234,1194,124) );//41
        sides.push_back( p0 + seg(1194,124,1144,42) );//42
        sides.push_back( p0 + seg(1144,42,1083,10) );//43
        sides.push_back( p0 + seg(1083,10,1012,7) );//44
        sides.push_back( p0 + seg(1012,7,921,61) );//45
        sides.push_back( p0 + seg(921,61,843,10) );//46
        sides.push_back( p0 + seg(843,10,760,3) );//47
        sides.push_back( p0 + seg(760,3,658,82) );//48
        cells.push_back(SidesToBoundaryMap(16));
        cells.back()[0] = {32,0};
        cells.back()[1] = {34,1};
        cells.back()[2] = {35,1};
        cells.back()[3] = {36,1};
        cells.back()[4] = {37,1};
        cells.back()[5] = {38,1};
        cells.back()[6] = {39,1};
        cells.back()[7] = {40,1};
        cells.back()[8] = {41,1};
        cells.back()[9] = {42,1};
        cells.back()[10] = {43,1};
        cells.back()[11] = {44,1};
        cells.back()[12] = {45,1};
        cells.back()[13] = {46,1};
        cells.back()[14] = {47,1};
        cells.back()[15] = {48,1};
        
//        cells.push_back(SidesToBoundaryMap(13));
//        cells.back()[0] = {24,0};
//        cells.back()[1] = {18,1};
        
        Real outer_dG = 10.*dG;
        outer_bnd = Boundary::get_ellipse({5e-2,2e-2},20e-2,20e-2,outer_dG,1.,0.);
        
        inner_dom_desc=cells;
        return;
    }   
//    if(n==4)
//    {
//        Real r1 = 1e-2;
//        Real r2 = 1.5e-2;
//        Real outer_dG = 10.*dG;
//        Boundary out_bnd(Boundary::get_circle({0.,0.},2.*r2,outer_dG,1.));
//        
//        vector<SidesToBoundaryMap> cells(0);
//        vector<Side> sides(0);
//        
//        sides.push_back(Segment({0.,-r1},{0.,-r2},dG,1.));
//        sides.push_back(Segment({0.,r2},{0.,r1},dG,1.));
//        sides.push_back(-Arc({0.,0.},r1,-M_PI/2.,M_PI/2.,dG,1.));
//        sides.push_back(Arc({0.,0.},r2,-M_PI/2.,M_PI/2.,dG,1.));
//        sides.push_back(-Arc({0.,0.},r1,M_PI/2.,3.*M_PI/2.,dG,1.));
//        sides.push_back(Arc({0.,0.},r2,M_PI/2.,3.*M_PI/2.,dG,1.));
//        
//        cells.push_back(SidesToBoundaryMap(4));
//        cells.push_back(SidesToBoundaryMap(4));
//        
//        cells[0][0] = {0,1};
//        cells[0][1] = {3,1};
//        cells[0][2] = {1,1};
//        cells[0][3] = {2,1};
//        
//        cells[1][0] = {0,0};
//        cells[1][1] = {4,1};
//        cells[1][2] = {1,0};
//        cells[1][3] = {5,1};
//        
//        MultiDomain mdom(out_bnd,sides,cells);
//        return mdom;
//    }
    
    cerr<<"MultiDomain n="<<n<<" unknown. I return a simple one."<<endl;
    get_multidomain_desc(1,param,outer_bnd,sides,inner_dom_desc);
    
}

void add_cells_cluster(const Parameters& param, Point& p, 
        vector<SidesToBoundaryMap>& cells, vector<Side>& sides)
{
    unsigned nx = param.P20_nx;
    unsigned ny = param.P20_ny;
    Real cl = param.P20_cell_l;
    Real cw = param.P20_cell_w;
    Real dG = param.P20_dG;
    bool vsmoothwave = param.P20_vsmoothwave;
    bool hsmoothwave = param.P20_hsmoothwave;
    Real vertamp = param.P20_vertamp;
    unsigned vertfreq = param.P20_vertfreq;
    Real horamp = param.P20_horamp;
    Real horlen = param.P20_horlen;
    Real horpos = param.P20_horpos;
    int horalt = param.P20_horalt;
    Real kappa_l = 1./param.P20_Rl;
    Real kappa_t = 1./param.P20_Rt;
    Real kappa_t_dir,kappa_t_perp;
    kappa_t_perp = kappa_t;
    if(param.P20_onlyperp)
        kappa_t_dir = 0.;
    else
        kappa_t_dir = kappa_t;
    
    unsigned cell_n = cells.size();
    unsigned sides_n = sides.size();
    
    cells.resize(cell_n + nx*ny);
    sides.resize(sides_n + 4*nx*ny-nx*(ny-1)-ny*(nx-1));
    
    // cells have length cl
    // cells have height cw
     
    Side v11;
    if(vsmoothwave || vertamp==0. || vertfreq==0)
        v11 = Wave1({0,0},{0.,cw},vertfreq,vertamp,dG,1.,kappa_l);
    else
        v11 = SquaredWave1({0,0},{0.,cw},vertfreq,vertamp,dG,1.,kappa_l,kappa_l); 
     
    for(unsigned j=0;j<ny;j++)
    {
        sides[sides_n+(nx+1)*j] = p+Point(0.,j*cw)+(-v11);
        for(unsigned i=1;i<=nx;i++)
                sides[sides_n+(nx+1)*j+i] = p+Point(i*cl,j*cw)+v11;
    }
    
    if(horalt<=1) 
    {
        Side h11,h12;
        if(horalt==-1) // a flat side (segment) with constant permeability
        {
            h11 = Segment({0,0},{cl,0.},dG,1.,kappa_t);
            h12 = h11;
        }
        else if(horalt>=0)
        {
            //a Segment h11a followed by a Side h11b and another Segment h11c
            //h11b can be a Segment, a Wave or a nonsmooth wave (SquaredWave)
            //permeability will be non zero only in h11b
            //h11b is centered in horpors*cl
            Side h11a,h11b,h11c;
            if(horamp>0.)
            {
                if(hsmoothwave)
                    h11b = Wave2({horpos*cl-horlen/2.,0.},{horpos*cl+horlen/2.,0},1,horamp,dG,1.,kappa_t);
                else
                    h11b = SquaredWave2({horpos*cl-horlen/2.,0.},{horpos*cl+horlen/2.,0},1,horamp,dG,1.,kappa_t_dir,kappa_t_perp);
            }
            else
                h11b = Segment({horpos*cl-horlen/2.,0},{horpos*cl+horlen/2.,0.},dG,1.,kappa_t);
            h11a = Segment({0.,0.},{horpos*cl-horlen/2.,0},dG,1.,0.);
            h11c = Segment({horpos*cl+horlen/2.,0.},{cl,0.,},dG,1.,0.);        
            h11 = Side({h11a,h11b,h11c});

            //if horalt==1 we alternate the position of h11b to force the flux to do
            //a zigzag
            if(horalt==0)
                h12 = h11;
            else
            {
                Side h12a,h12b,h12c;
                horpos = 1.-horpos;
                if(horamp>0.)
                    if(hsmoothwave)
                        h12b = Wave2({horpos*cl-horlen/2.,0.},{horpos*cl+horlen/2.,0},1,horamp,dG,1.,kappa_t);
                    else
                        h12b = SquaredWave2({horpos*cl-horlen/2.,0.},{horpos*cl+horlen/2.,0},1,horamp,dG,1.,kappa_t_dir,kappa_t_perp);
                else
                    h12b = Segment({horpos*cl-horlen/2.,0},{horpos*cl+horlen/2.,0.},dG,1.,kappa_t);
                h12a = Segment({0.,0.},{horpos*cl-horlen/2.,0},dG,1.,0.);
                h12c = Segment({horpos*cl+horlen/2.,0.},{cl,0.,},dG,1.,0.);
                h12 = Side({h12a,h12b,h12c});
                horpos = 1.-horpos;
            }
        }

        for(unsigned j=0;j<ny;j++)
        for(unsigned i=0;i<nx;i++)
        {
            if(j%2==0)
                sides[sides_n+ny*(nx+1)+nx*j+i] = p+Point(i*cl,j*cw)+h11;
            else
                sides[sides_n+ny*(nx+1)+nx*j+i] = p+Point(i*cl,j*cw)+h12;
        }

        if(ny%2==0)
            for(unsigned i=0;i<nx;i++)
                sides[sides_n+ny*(nx+1)+nx*ny+i] = p+Point(i*cl,ny*cw)+(-h11);
        else
            for(unsigned i=0;i<nx;i++)
                sides[sides_n+ny*(nx+1)+nx*ny+i] = p+Point(i*cl,ny*cw)+(-h12);
    }
    else if(horalt==2)
    {
        unsigned myseed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::minstd_rand0 generator1;
        generator1.seed(myseed);
        std::uniform_real_distribution<double> distribution1(horlen/cl,1.-horlen/cl);
        
        std::uniform_real_distribution<double> distribution2(0.,1.);
        std::minstd_rand0 generator2;
        myseed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        generator2.seed(myseed);

        Real horprob = param.P20_horprob;
        Real kappa_mult;
        
        Side h11;        
        Side h11a,h11b,h11c;        
        for(unsigned j=0;j<ny;j++)
        for(unsigned i=0;i<nx;i++)
        {
            horpos = distribution1(generator1);
            if(distribution2(generator2)<=horprob)
                kappa_mult=0.;
            else
                kappa_mult=1.;
            
            if(horamp>0.)
            {
                if(hsmoothwave)
                    h11b = Wave2({horpos*cl-horlen/2.,0.},{horpos*cl+horlen/2.,0},1,horamp,dG,1.,kappa_mult*kappa_t);
                else
                    h11b = SquaredWave2({horpos*cl-horlen/2.,0.},{horpos*cl+horlen/2.,0},1,horamp,dG,1.,kappa_mult*kappa_t_dir,kappa_mult*kappa_t_perp);
            }
            else
                h11b = Segment({horpos*cl-horlen/2.,0},{horpos*cl+horlen/2.,0.},dG,1.,kappa_mult*kappa_t);
            h11a = Segment({0.,0.},{horpos*cl-horlen/2.,0},dG,1.,0.);
            h11c = Segment({horpos*cl+horlen/2.,0.},{cl,0.,},dG,1.,0.);        
            h11 = Side({h11a,h11b,h11c});
            sides[sides_n+ny*(nx+1)+nx*j+i] = p+Point(i*cl,j*cw)+h11;
        }

        for(unsigned i=0;i<nx;i++)
        {
            horpos = distribution1(generator1);
            if(horamp>0.)
            {
                if(hsmoothwave)
                    h11b = Wave2({horpos*cl-horlen/2.,0.},{horpos*cl+horlen/2.,0},1,horamp,dG,1.,kappa_t);
                else
                    h11b = SquaredWave2({horpos*cl-horlen/2.,0.},{horpos*cl+horlen/2.,0},1,horamp,dG,1.,kappa_t_dir,kappa_t_perp);
            }
            else
                h11b = Segment({horpos*cl-horlen/2.,0},{horpos*cl+horlen/2.,0.},dG,1.,kappa_t);
            h11a = Segment({0.,0.},{horpos*cl-horlen/2.,0},dG,1.,0.);
            h11c = Segment({horpos*cl+horlen/2.,0.},{cl,0.,},dG,1.,0.);        
            h11 = Side({h11a,h11b,h11c});
            sides[sides_n+ny*(nx+1)+nx*ny+i] = p+Point(i*cl,ny*cw)+(-h11);
        }
    }
    
    
    for(unsigned j=0;j<ny;j++)
    for(unsigned i=0;i<nx;i++)
    {
        if(i>0 && j<ny-1)
        {
            cells[cell_n+nx*j+i] = SidesToBoundaryMap(4);
            cells[cell_n+nx*j+i][0] = {sides_n+ j*(nx+1)+i,0};
            cells[cell_n+nx*j+i][1] = {sides_n+ ny*(nx+1)+j*nx+i,1};
            cells[cell_n+nx*j+i][2] = {sides_n+ j*(nx+1)+i+1,1};
            cells[cell_n+nx*j+i][3] = {sides_n+ ny*(nx+1)+(j+1)*nx+i,0};
        }        
        else if(i==0 && j<ny-1)
        {
            cells[cell_n+nx*j+i] = SidesToBoundaryMap(4);
            cells[cell_n+nx*j+i][0] = {sides_n+ j*(nx+1)+i,1};
            cells[cell_n+nx*j+i][1] = {sides_n+ ny*(nx+1)+j*nx+i,1};
            cells[cell_n+nx*j+i][2] = {sides_n+ j*(nx+1)+i+1,1};
            cells[cell_n+nx*j+i][3] = {sides_n+ ny*(nx+1)+(j+1)*nx+i,0};
        }
        else if(i>0 && j==ny-1)
        {
            cells[cell_n+nx*j+i] = SidesToBoundaryMap(4);
            cells[cell_n+nx*j+i][0] = {sides_n+ j*(nx+1)+i,0};
            cells[cell_n+nx*j+i][1] = {sides_n+ ny*(nx+1)+j*nx+i,1};
            cells[cell_n+nx*j+i][2] = {sides_n+ j*(nx+1)+i+1,1};
            cells[cell_n+nx*j+i][3] = {sides_n+ ny*(nx+1)+(j+1)*nx+i,1};
        }
        else if(i==0 && j==ny-1)
        {
            cells[cell_n+nx*j+i] = SidesToBoundaryMap(4);
            cells[cell_n+nx*j+i][0] = {sides_n+ j*(nx+1)+i,1};
            cells[cell_n+nx*j+i][1] = {sides_n+ ny*(nx+1)+j*nx+i,1};
            cells[cell_n+nx*j+i][2] = {sides_n+ j*(nx+1)+i+1,1};
            cells[cell_n+nx*j+i][3] = {sides_n+ ny*(nx+1)+(j+1)*nx+i,1};
        }
    }
}

//void add_cells_cluster_with_holes(const Real cl, const Real cw, 
//        const Real pl, const Real pw, Point& p, 
//        unsigned nx, unsigned ny, const Real dGl, const Real dGw,
//        vector<SidesToBoundaryMap>& cells, vector<Side>& sides)
//{
//    unsigned cell_n = cells.size();
//    unsigned sides_n = sides.size();
//    
//    const unsigned n_corner_sides = 12;
//    
//    cells.resize(cell_n + nx*ny);
//    sides.resize(sides_n + 4*nx*ny-nx*(ny-1)-ny*(nx-1)
//                 +n_corner_sides*nx*ny);
//    
//    // cells have length cl
//    // cells have height cw
//    // horizontal gap junctions have length gkl
//    // vertical gap junctions have length gkw
//    const Real gjl = pl*cl;
//    const Real gjw = pw*cw;
//    // h1,h2,v1,v2 are gap junctions
//    Side h1 = Point((cl-gjl)/2.,0.)+Segment({0,0},{gjl,0.},dGl,1.);
//    Side h2 = Point(-(cl-gjl)/2.,0.)+Segment({0,0},{-gjl,0.},dGl,1.);
//    Side v1 = Point(0.,(cw-gjw)/2.)+Segment({0,0},{0.,gjw},dGw,1.);
//    Side v2 = Point(0.,-(cw-gjw)/2.)+Segment({0,0},{0.,-gjw},dGw,1.);
//    
//    //first, we add the gap junctions
//    for(unsigned j=0;j<ny;j++)
//    {
//        sides[sides_n+(nx+1)*j] = p+Point(0.,(j+1)*cw)+v2;
//        for(unsigned i=1;i<=nx;i++)
//            sides[sides_n+(nx+1)*j+i] = p+Point(i*cl,j*cw)+v1;
//    }
//    
//    for(unsigned j=0;j<ny;j++)
//    for(unsigned i=0;i<nx;i++)
//        sides[sides_n+ny*(nx+1)+nx*j+i] = p+Point(i*cl,j*cw)+h1;
//    
//    for(unsigned i=0;i<nx;i++)
//        sides[sides_n+ny*(nx+1)+nx*ny+i] = p+Point((i+1)*cl,ny*cw)+h2;
//    
//    //now we add the transmembrane sides
//    Real reml=(cl-gjl)/2.;
//    Real remw=(cw-gjw)/2.;
//    vector<Side> corners;
//    if(n_corner_sides==12)
//    {
//        Real dl = 0.05;
//        Real dw = 0.05;        
//        Real dgjl = dl*cw;
//        Real dgjw = dw*cl;
//        Real r1 = reml-dgjw;
//        Real r2 = remw-dgjl;
//        if(r1<=0. || r2<=0.)
//        {
//            cout<<"ERROR constructing the cells: dl or dw is too large"<<endl;
//            return;
//        }
//
//        corners.push_back( Segment({0.,remw},{dgjw,remw},dGl,1.) );
//        corners.push_back( EllipseArc({reml,remw},r1,r2,M_PI,3.*M_PI/2.,dGl,1.) );
//        corners.push_back( Segment({reml,dgjl},{reml,0.},dGl,1.) );
//
//        corners.push_back( Segment({cl-reml,0.},{cl-reml,dgjl},dGl,1.) );
//        corners.push_back( EllipseArc({cl-reml,remw},r1,r2,-M_PI/2.,0.,dGl,1.) );
//        corners.push_back( Segment({cl-dgjw,remw},{cl,remw},dGl,1.) );
//
//        corners.push_back( Segment({cl,cw-remw},{cl-dgjw,cw-remw},dGl,1.) );
//        corners.push_back( EllipseArc({cl-reml,cw-remw},r1,r2,0.,M_PI/2.,dGl,1.) );
//        corners.push_back( Segment({cl-reml,cw-dgjl},{cl-reml,cw},dGl,1.) );
//
//        corners.push_back( Segment({reml,cw},{reml,cw-dgjl},dGl,1.) );
//        corners.push_back( EllipseArc({reml,cw-remw},r1,r2,M_PI/2.,M_PI,dGl,1.) );
//        corners.push_back( Segment({dgjw,cw-remw},{0.,cw-remw},dGl,1.) );
//    }
//    else if(n_corner_sides==4)
//    {               
//        corners.push_back( -EllipseArc({0.,0.},reml,remw,0.,M_PI/2.,dGl,1.) );
//        corners.push_back( -EllipseArc({cl,0.},reml,remw,M_PI/2.,M_PI,dGl,1.) );
//        corners.push_back( -EllipseArc({cl,cw},reml,remw,M_PI,3.*M_PI/2.,dGl,1.) );
//        corners.push_back( -EllipseArc({0.,cw},reml,remw,3.*M_PI/2.,2.*M_PI,dGl,1.) );
//    }
//    else
//    {
//        cout<<"Error: n_corner_sides must be either 4 or 12"<<endl;
//        return;
//    }
//    
//    
//    unsigned start = sides_n+ny*(nx+1)+nx*(ny+1);
//    for(unsigned j=0;j<ny;j++)
//    for(unsigned i=0;i<nx;i++)
//    for(unsigned k=0;k<n_corner_sides;k++)
//        sides[start+n_corner_sides*(j*nx+i)+k] = p+Point(i*cl,j*cw)+corners[k];
//        
//    
//    for(unsigned j=0;j<ny;j++)
//    for(unsigned i=0;i<nx;i++)
//    {
//        cells[cell_n+nx*j+i] = SidesToBoundaryMap(4+n_corner_sides);
//        if(i>0 && j<ny-1)
//        {
//            
//            cells[cell_n+nx*j+i][0] = {sides_n+ j*(nx+1)+i,0};
//            cells[cell_n+nx*j+i][1] = {sides_n+ ny*(nx+1)+j*nx+i,1};
//            cells[cell_n+nx*j+i][2] = {sides_n+ j*(nx+1)+i+1,1};
//            cells[cell_n+nx*j+i][3] = {sides_n+ ny*(nx+1)+(j+1)*nx+i,0};            
//        }        
//        else if(i==0 && j<ny-1)
//        {
//            cells[cell_n+nx*j+i][0] = {sides_n+ j*(nx+1)+i,1};
//            cells[cell_n+nx*j+i][1] = {sides_n+ ny*(nx+1)+j*nx+i,1};
//            cells[cell_n+nx*j+i][2] = {sides_n+ j*(nx+1)+i+1,1};
//            cells[cell_n+nx*j+i][3] = {sides_n+ ny*(nx+1)+(j+1)*nx+i,0};
//        }
//        else if(i>0 && j==ny-1)
//        {
//            cells[cell_n+nx*j+i][0] = {sides_n+ j*(nx+1)+i,0};
//            cells[cell_n+nx*j+i][1] = {sides_n+ ny*(nx+1)+j*nx+i,1};
//            cells[cell_n+nx*j+i][2] = {sides_n+ j*(nx+1)+i+1,1};
//            cells[cell_n+nx*j+i][3] = {sides_n+ ny*(nx+1)+(j+1)*nx+i,1};
//        }
//        else if(i==0 && j==ny-1)
//        {
//            cells[cell_n+nx*j+i][0] = {sides_n+ j*(nx+1)+i,1};
//            cells[cell_n+nx*j+i][1] = {sides_n+ ny*(nx+1)+j*nx+i,1};
//            cells[cell_n+nx*j+i][2] = {sides_n+ j*(nx+1)+i+1,1};
//            cells[cell_n+nx*j+i][3] = {sides_n+ ny*(nx+1)+(j+1)*nx+i,1};
//        }
//        
//        for(unsigned k=0;k<n_corner_sides;k++)
//            cells[cell_n+nx*j+i][4+k] = {start+n_corner_sides*(j*nx+i)+k,1};
//    }
//}

//Side squared_piramid(const Point a, const Point b, 
//                     const Real H, const Real dG)
//{
//    Point dir = b-a;
//    Real ls = dir.norm()/6.;
//    Real hs = H/2.;
//    dir = dir/dir.norm();
//    Point perp(-dir(1),dir(0));
//    
//    Side h = Segment({0.,0.},ls*dir,dG,1.);
//    Side v = Segment({0.,0.},hs*perp,dG,1.);
//    Side vr = Segment({0.,0.},-hs*perp,dG,1.);
//    
//    vector<Side> sides(10);
//    sides[0] = a+h;
//    sides[1] = sides[0].end()+v;
//    sides[2] = sides[1].end()+h;
//    sides[3] = sides[2].end()+v;
//    sides[4] = sides[3].end()+h;
//    sides[5] = sides[4].end()+h;
//    sides[6] = sides[5].end()+vr;
//    sides[7] = sides[6].end()+h;
//    sides[8] = sides[7].end()+vr;
//    sides[9] = sides[8].end()+h;
//    
////    Real ls = L/3.;
////    Real hs = H;    
////    Side h = Segment({0.,0.},{ls,0.},dG);
////    Side v = Segment({0.,0.},{0.,hs},dG);
////    Side vr = Segment({0.,0.},{0.,-hs},dG);
////    
////    vector<Side> sides(5);
////    sides[0] = h;
////    sides[1] = sides[0].end()+v;
////    sides[2] = sides[1].end()+h;
////    sides[3] = sides[2].end()+vr;
////    sides[4] = sides[3].end()+h;
//    
////    Real ls = L/4.;
////    Real hs = H;    
////    Side h = Segment({0.,0.},{ls,0.},dG);
////    Side v = Segment({0.,0.},{0.,hs},dG);
////    Side vr = Segment({0.,0.},{0.,-hs},dG);
////    
////    vector<Side> sides(8);
////    sides[0] = h;
////    sides[1] = sides[0].end()+v;
////    sides[2] = sides[1].end()+h;
////    sides[3] = sides[2].end()+vr;
////    sides[4] = sides[3].end()+vr;
////    sides[5] = sides[4].end()+h;
////    sides[6] = sides[5].end()+v;
////    sides[7] = sides[6].end()+h;
//    
//    return Side(sides);
//}