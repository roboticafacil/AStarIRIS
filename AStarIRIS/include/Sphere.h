#ifndef SPHERE
#define SPHERE
#include "fusion.h"
#include "ConicSet.h"

using namespace mosek::fusion;
using namespace monty;


class Sphere :  public ConicSet
{
public:
    double radius;
    Sphere(const int& n);
    Sphere(const Sphere& sphere);
    Sphere(const Eigen::VectorXd& center,const double &radius);
    ~Sphere();
    Sphere& operator=(const Sphere& other);
    double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    double closestPoint(const Eigen::VectorXd& p_in);
    double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
    bool isInside(const Eigen::VectorXd& p, const double &tol=1.e-4);
    bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol=1.e-4);
    bool isInsideSeparatingHyperplanes(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const double& tol = 1.e-4);
    Eigen::VectorXd getCentroid();
    Eigen::MatrixXd getBoundingBox();
    void allocateClosestPointEllipsoidSolver();
private:
    std::shared_ptr<ndarray<double, 1>> c_ptr;
    Model::t MExpandingEllipsoid;
    Variable::t xExpandingEllipsoid;
    Variable::t alphaExpandingEllipsoid;
    bool solverExpandingEllipsoidAllocated = false;
};
#endif