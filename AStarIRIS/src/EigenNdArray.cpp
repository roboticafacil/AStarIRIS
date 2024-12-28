#include<Eigen/Dense>
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

std::shared_ptr<ndarray<double, 2>> Eigen2NdArray(const Eigen::MatrixXd& A)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Arow(A);
    std::shared_ptr<ndarray<double, 2>> A_ptr(new ndarray<double, 2>(Arow.data(), shape(A.rows(), A.cols())));
    return A_ptr;
}

std::shared_ptr<ndarray<double, 1>> Eigen2NdArray(const Eigen::VectorXd& b)
{
    std::shared_ptr<ndarray<double, 1>> b_ptr(new ndarray<double, 1>(b.data(), shape(b.rows())));
    return b_ptr;
}