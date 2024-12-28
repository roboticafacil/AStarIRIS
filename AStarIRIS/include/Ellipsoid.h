#ifndef ELLIPSOID
#define ELLIPSOID
#include "fusion.h"
#include "ConicSet.h"

using namespace mosek::fusion;
using namespace monty;

class Ellipsoid: public ConicSet
{
public:
    Eigen::MatrixXd C;
    Ellipsoid(const int& n);
    Ellipsoid(const Eigen::MatrixXd& C, const Eigen::VectorXd& d);
    Ellipsoid(const Ellipsoid& ellipsoid);
    ~Ellipsoid();
    Ellipsoid& operator=(const Ellipsoid& other);
    void print();
    double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    double closestPoint(const Eigen::VectorXd& p_in);
    double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
    bool isInside(const Eigen::VectorXd& p);
    bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi);
    Eigen::VectorXd getCentroid();
    Eigen::MatrixXd getCInv();
    void allocateSolver();
    Eigen::MatrixXd getBoundingBox();
    void tangentHyperplane(const Eigen::VectorXd& p_in, Eigen::VectorXd& ai, double& bi);
    void allocateIsInsideSeparatingHyperplaneSolver();
    void allocateClosestPointEllipsoidSolver();
private:
    std::shared_ptr<ndarray<double, 2>> Cinv_ptr;
    std::shared_ptr<ndarray<double, 1>> centroid_ptr;
    Model::t M;
    Variable::t x;
    Variable::t alpha;
    bool solverAllocated=false;
    Eigen::MatrixXd Cinv;
    bool inverseComputed = false;
    Model::t MIsInsideSeparatingHyperplane;
    Variable::t xIsInsideSeparatingHyperplane;
    bool solverIsInsideSeparatingHyperplaneAllocated = false;
    Model::t MExpandingEllipsoid;
    Variable::t xExpandingEllipsoid;
    Variable::t alphaExpandingEllipsoid;
    bool solverExpandingEllipsoidAllocated = false;
};
#endif