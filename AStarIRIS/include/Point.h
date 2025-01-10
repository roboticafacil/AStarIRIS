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
	std::ostream& print(std::ostream& out, const std::string& pointName);
	double distance(const Eigen::VectorXd& point);
	double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
	double closestPoint(const Eigen::VectorXd& p_in);
	double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
	bool isInside(const Eigen::VectorXd& p, const double &tol= 1.e-4);
	bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol=1.e-4);
	Eigen::VectorXd getCentroid();
	Eigen::MatrixXd getBoundingBox();
};
#endif