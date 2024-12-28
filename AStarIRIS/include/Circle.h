#ifndef CIRCLE
#define CIRCLE
#include "fusion.h"
#include "ConicSet.h"

using namespace mosek::fusion;
using namespace monty;


class Circle :  public ConicSet
{
public:
    double radius;
	Circle(const int& n);
    Circle(const Circle& circle);
	Circle(const Eigen::VectorXd& center,const double &radius);
    ~Circle();
    Circle& operator=(const Circle& other);
    double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    double closestPoint(const Eigen::VectorXd& p_in);
    double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
    bool isInside(const Eigen::VectorXd& p);
    bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi);
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