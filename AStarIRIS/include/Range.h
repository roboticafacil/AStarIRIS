#ifndef RANGE
#define RANGE
#include <Eigen/Dense>
#include "fusion.h"

using namespace monty;

class Range
{
public:
    void setRange(const int& n,const double& lb, const double& ub);
    std::pair<double,double> getBounds();
    std::pair<std::shared_ptr<ndarray<double, 2>>,std::shared_ptr<ndarray<double, 1>>> getConstraints();
    void getConstraints(Eigen::MatrixXd& A, Eigen::VectorXd &b);
protected:
    Range();
    static Range* _range;
    double _lb;
    double _ub;
    std::shared_ptr<ndarray<double, 2>> A_ptr;
    std::shared_ptr<ndarray<double, 1>> b_ptr;
public:
    Range(Range& other) = delete;
    void operator=(const Range&) = delete;
    static Range* getInstance();
private:
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
};
#endif
