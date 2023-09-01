#ifndef MAINHEADER_H
#define MAINHEADER_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <deque>
#include <random>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <stdio.h>
#include <string.h>

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

using namespace std;

typedef double Real;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::ArrayXd Array;
typedef Eigen::SparseMatrix<double> SpMatrix;

//#include <autodiff/forward/real.hpp>
//#include <autodiff/forward/real/eigen.hpp>
//typedef autodiff::VectorXreal ADVector;
//typedef autodiff::ArrayXreal ADArray;
//typedef autodiff::real ADReal;

enum Equation {ODE,D_SDE,JD_SDE};
enum Controller {NONE=0,I=1,PI=2,PPI=3}; 
// Integral, Proportional Integral, Predictive Proportional Integral

enum IonicModel {NO_GATING_VARS, RogersMcCulloch, 
            MitchellSchaeffer, BuenoCherryFenton};


#endif /* MAINHEADER_H */

