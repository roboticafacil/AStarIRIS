#ifndef CONIC_SET
#define CONIC_SET
#include<Eigen/Dense>

typedef enum ConicType
{
    point = 0,
    circle = 1,
    ellipsoid = 2,
    polyhedron = 3,
    obstacleCircularRobot = 4
} ConicType;

class Ellipsoid;

class ConicSet
{
public:
    const int n;
protected:
        bool centroidComputed =false;
        Eigen::VectorXd centroid;
        bool bbComputed = false;
        Eigen::MatrixXd bb;
public:
    ConicSet(const int& n);
    ConicSet(const ConicSet& conicSet);
    ~ConicSet();
    virtual double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out) = 0;
    virtual double closestPoint(const Eigen::VectorXd& p_in)=0;
    virtual double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out) = 0;
    virtual bool isInside(const Eigen::VectorXd& p, const double &tol=0.)=0;
    virtual bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol=0.) = 0;
    virtual Eigen::VectorXd getCentroid() = 0;
    virtual Eigen::MatrixXd getBoundingBox() = 0;
};
#endif