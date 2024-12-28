#ifndef POINT
#define POINT
#include<Eigen/Dense>
#include "ConicSet.h"

class Point : public ConicSet
{
public:
	Eigen::VectorXd p;
	Point(const Eigen::VectorXd& p);
	Point(const Point& point);
	void print();
	double distance(const Eigen::VectorXd& point);
	double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
	double closestPoint(const Eigen::VectorXd& p_in);
	double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
	bool isInside(const Eigen::VectorXd& p);
	bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi);
	Eigen::VectorXd getCentroid();
	Eigen::MatrixXd getBoundingBox();
};
#endif