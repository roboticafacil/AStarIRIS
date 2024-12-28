#ifndef EIGEN_ND_ARRAY
#define EIGEN_ND_ARRAY
#include<Eigen/Dense>
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

std::shared_ptr<ndarray<double, 2>> Eigen2NdArray(const Eigen::MatrixXd& A);
std::shared_ptr<ndarray<double, 1>> Eigen2NdArray(const Eigen::VectorXd& b);

#endif