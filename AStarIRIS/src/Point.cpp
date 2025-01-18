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

std::ostream& Point::print(std::ostream& out, const std::string& pointName)
{
	out << pointName << "=[";
	for (int i = 0; i < this->p.rows(); i++)
	{
		if (i < (this->p.rows() - 1))
			out << this->p(i) << ";";
		else
			out << this->p(i);
	}
	out << "];" << std::endl;
	return out;
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

bool Point::isInside(const Eigen::VectorXd& p, const double &tol)
{
	return (this->p-p).norm()<=tol;
}

bool Point::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol)
{
	return ((ai.dot(p) - bi*ai.norm())<=tol);
}

bool Point::isInsideSeparatingHyperplanes(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const double& tol)
{
	std::cout << "Point: isInsideSeparatingHyperplanes not implemented yet!!" << std::endl;
	//Hay que comprobar esto...
	return ((((A * this->p) - b).array()) <= tol).all();
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
