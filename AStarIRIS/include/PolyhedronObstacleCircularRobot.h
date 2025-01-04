#ifndef POLYHEDRON_OBSTACLE_CIRCULAR_ROBOT
#define POLYHEDRON_OBSTACLE_CIRCULAR_ROBOT
#include<Eigen/Dense>
#include "PolyhedronV.h"
#include "ConicSet.h"

class PolyhedronObstacleCircularRobot : public PolyhedronV//, public ConicSet
{
public:
    double r;
public:
    PolyhedronObstacleCircularRobot(const Eigen::MatrixXd& v, const double& r);
    PolyhedronObstacleCircularRobot(const PolyhedronObstacleCircularRobot& polyhedronCircularRobot);
    double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    double closestPoint(const Eigen::VectorXd& p_in);
    bool isInside(const Eigen::VectorXd& p, const double &tol= 1.e-4);
    bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol= 1.e-4);
    Ellipsoid inscribedEllipsoid();
    Circle inscribedCircle();
    Eigen::VectorXd getCentroid();
    double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
    void allocateClosestPointSolver();
    void allocateClosestPointEllipsoidSolver();
    void allocateIsInsideSeparatingHyperplaneSolver();
    Eigen::MatrixXd getBoundingBox();
private:
    Variable::t chi1;
    Variable::t chi2;
    Variable::t chi1ClosestPointEllipsoid;
    Variable::t chi2ClosestPointEllipsoid;
    Variable::t chi1IsInsideSeparatingHyperplane;
    Variable::t chi2IsInsideSeparatingHyperplane;
    //Parameter::t radius;
};
#endif