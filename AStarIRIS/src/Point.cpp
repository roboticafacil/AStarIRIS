#include<Eigen/Dense>
#include "Point.h"
#include "Ellipsoid.h"
#include <iostream>

Point::Point(const Eigen::VectorXd& p): ConicSet(p.rows()) //ConicSet(p.rows(),false)
{
	this->p = p;
}

Point::Point(const Point& point): ConicSet(point), p(point.p) 
{
}

void Point::print()
{
	std::cout << "p=" << this->p << std::endl;
}

double Point::distance(const Eigen::VectorXd& point)
{
	return (p - point).norm();
}

double Point::closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out)
{
	p_out = this->p;
	return this->distance(p_in);
}
double Point::closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out)
{
	Eigen::VectorXd pg=ellipsoid.getCInv()* (this->p - ellipsoid.getCentroid());
	return pg.norm();
}
double Point::closestPoint(const Eigen::VectorXd& p_in)
{
	return this->distance(p_in);
}

bool Point::isInside(const Eigen::VectorXd& p)
{
	return (this->p-p).norm()<1e-7;
}

bool Point::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi)
{
	return ((ai.dot(p) - bi*ai.norm())<=0.);
}

Eigen::VectorXd Point::getCentroid()
{
	if (!centroidComputed)
	{
		this->centroidComputed = true;
		this->centroid = this->p;
	}
	return this->centroid;
}

Eigen::MatrixXd Point::getBoundingBox()
{
	if (!this->bbComputed)
	{
		Eigen::MatrixXd bb(this->n, 2);
		bb.col(0) = this->p;
		bb.col(1) = this->p;
		this->bbComputed = true;
	}
	return bb;
}
