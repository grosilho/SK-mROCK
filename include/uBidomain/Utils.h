#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <iostream>
#include <fstream>

namespace Eigen {
// taken from https://stackoverflow.com/a/25389481/11927397
template<class Matrix>
inline void write_binary(std::ofstream& out, const Matrix& matrix){
//    std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    if(out.is_open()) 
    {
        typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
        out.write(reinterpret_cast<const char*>(&rows), sizeof(typename Matrix::Index));
        out.write(reinterpret_cast<const char*>(&cols), sizeof(typename Matrix::Index));
        out.write(reinterpret_cast<const char*>(matrix.data()), rows*cols*static_cast<typename Matrix::Index>(sizeof(typename Matrix::Scalar)) );
//        out.close();
    }
    else 
        std::cout << "Can not write to file."<< std::endl;
    
}

template<class Matrix>
inline void read_binary(std::ifstream& in, Matrix& matrix){
//    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (in.is_open()) {
        typename Matrix::Index rows=0, cols=0;
        in.read(reinterpret_cast<char*>(&rows),sizeof(typename Matrix::Index));
        in.read(reinterpret_cast<char*>(&cols),sizeof(typename Matrix::Index));
        matrix.resize(rows, cols);
        in.read(reinterpret_cast<char*>(matrix.data()), rows*cols*static_cast<typename Matrix::Index>(sizeof(typename Matrix::Scalar)) );
//        in.close();
    }
    else 
        std::cout << "Can not open binary matrix file."<< std::endl;
    
}

}

template<typename scalar>
inline void write_binary(std::ofstream& out, const vector<scalar>& vec)
{
    if(out.is_open()) 
    {
        unsigned N = vec.size();
        out.write(reinterpret_cast<const char*>(&N), sizeof(unsigned));       
        out.write(reinterpret_cast<const char*>(vec.data()), N*sizeof(scalar));
    }
    else 
        std::cout << "Can not write to file."<< std::endl;
    
}

template<typename scalar>
inline void read_binary(std::ifstream& in, vector<scalar>& vec)
{
    if(in.is_open()) 
    {
        unsigned N;
        in.read(reinterpret_cast<char*>(&N), sizeof(unsigned));
        vec.clear();
        vec.resize(N);        
        in.read(reinterpret_cast<char*>(vec.data()), N*sizeof(scalar));        
    }
    else 
        std::cout << "Can not read from file."<< std::endl;
    
}

template<typename scalar, unsigned N>
inline void write_binary(std::ofstream& out, const vector<array<scalar,N>>& vec_arr)
{
    if(out.is_open()) 
    {
        unsigned M = vec_arr.size();
        unsigned Nb = vec_arr[0].size();
        out.write(reinterpret_cast<const char*>(&M), sizeof(unsigned));
        out.write(reinterpret_cast<const char*>(&Nb), sizeof(unsigned));
        for(auto& arr: vec_arr)
            out.write(reinterpret_cast<const char*>(arr.data()), N*sizeof(scalar));
    }
    else 
        std::cout << "Can not write to file."<< std::endl;
    
}

template<typename scalar, unsigned N>
inline void read_binary(std::ifstream& in, vector<array<scalar,N>>& vec_arr)
{
    if(in.is_open()) 
    {
        unsigned M, N_check;
        in.read(reinterpret_cast<char*>(&M), sizeof(unsigned));
        in.read(reinterpret_cast<char*>(&N_check), sizeof(unsigned));
        if(N!=N_check)
        {
            std::cout<<"ERROR: N not the same as in input file."<<endl;
            return;
        }
        
        vec_arr.clear();
        vec_arr.resize(M);
        for(unsigned i=0;i<M;i++)
            in.read(reinterpret_cast<char*>(vec_arr[i].data()), N*sizeof(scalar));        
    }
    else 
        std::cout << "Can not read from file."<< std::endl;
    
}

template<typename scalar>
inline void write_scalar_binary(std::ofstream& out, const scalar& x)
{
    if(out.is_open())     
        out.write(reinterpret_cast<const char*>(&x), sizeof(scalar));            
    else 
        std::cout << "Can not write to file."<< std::endl;    
}

template<typename scalar>
inline void read_scalar_binary(std::ifstream& in, scalar& x)
{
    if(in.is_open())     
        in.read(reinterpret_cast<char*>(&x), sizeof(scalar));            
    else 
        std::cout << "Can not read from file."<< std::endl;    
}

inline void write_binary(std::ofstream& out, const vector<tuple<unsigned,unsigned,unsigned,Real,Real>>& vec_tup)
{
    if(out.is_open()) 
    {
        unsigned N = vec_tup.size();
        out.write(reinterpret_cast<const char*>(&N), sizeof(unsigned));
        for(auto& tup: vec_tup)
        {
            write_scalar_binary(out,std::get<0>(tup));
            write_scalar_binary(out,std::get<1>(tup));
            write_scalar_binary(out,std::get<2>(tup));
            write_scalar_binary(out,std::get<3>(tup));
            write_scalar_binary(out,std::get<4>(tup));
        }
    }
    else 
        std::cout << "Can not write to file."<< std::endl;
    
}

inline void read_binary(std::ifstream& in, vector<tuple<unsigned,unsigned,unsigned,Real,Real>>& vec_tup)
{
    if(in.is_open()) 
    {
        unsigned N;
        in.read(reinterpret_cast<char*>(&N), sizeof(unsigned));        
        vec_tup.clear();
        vec_tup.resize(N);
        for(unsigned i=0;i<N;i++)
        {
            unsigned a,b,c;
            Real d,e;
            read_scalar_binary(in,a);
            read_scalar_binary(in,b);
            read_scalar_binary(in,c);
            read_scalar_binary(in,d);
            read_scalar_binary(in,e);
            vec_tup[i] = {a,b,c,d,e};
        }       
    }
    else 
        std::cout << "Can not read from file."<< std::endl;
    
}

#endif /* UTILS_H */

