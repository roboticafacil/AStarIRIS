#ifndef EIGEN_UTILS
#define EIGEN_UTITLS
#include<Eigen/Dense>
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

std::shared_ptr<ndarray<double, 2>> Eigen2NdArray(const Eigen::MatrixXd& A);
std::shared_ptr<ndarray<double, 1>> Eigen2NdArray(const Eigen::VectorXd& b);
double nchoosek(const int n, const int k);
void nchoosek(const Eigen::VectorXi& V, const int k, Eigen::MatrixXi& U);
void unique(Eigen::MatrixXd& A);
bool areRowsEqual(const Eigen::RowVectorXd& row1, const Eigen::RowVectorXd& row2, const double& tolerance);
Eigen::MatrixXd removeDuplicateRows(const Eigen::MatrixXd& matrix, const double& tolerance);

#endif