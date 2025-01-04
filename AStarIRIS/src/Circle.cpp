#include "Circle.h"
#include "Ellipsoid.h"
#include "ConicSet.h"
#include "Range.h"
#include "EigenNdArray.h"


Circle::Circle(const int& n): radius(0.), ConicSet(n) //ConicSet(n,false)
{
}

Circle::Circle(const Circle& circle): radius(circle.radius), ConicSet(circle)
{
}

Circle::Circle(const Eigen::VectorXd& center, const double& radius): radius(radius), ConicSet(center.rows()) //ConicSet(center.rows(),false)
{
	this->centroid = center;
}

Circle::~Circle()
{
	//auto _M = finally([&]() { M->dispose(); });
}

Circle& Circle::operator=(const Circle& other) {
	this->radius = other.radius;
	this->alphaExpandingEllipsoid = other.alphaExpandingEllipsoid;
	this->bb = other.bb;
	this->bbComputed = other.bbComputed;
	this->centroid = other.centroid;
	this->centroidComputed = other.centroidComputed;
	this->MExpandingEllipsoid = other.MExpandingEllipsoid;
	this->solverExpandingEllipsoidAllocated = other.solverExpandingEllipsoidAllocated;
	this->xExpandingEllipsoid = other.xExpandingEllipsoid;
	this->c_ptr = other.c_ptr;
	return *this;
}

double Circle::closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out)
{
	std::shared_ptr<ndarray<double, 2>> Cinv_ptr = Eigen2NdArray(ellipsoid.getCInv());
	std::shared_ptr<ndarray<double, 1>> centroid_ptr = Eigen2NdArray(ellipsoid.getCentroid());
	this->allocateClosestPointEllipsoidSolver();
	if (!this->MExpandingEllipsoid->hasConstraint("ellipsoid"))
	{
		this->MExpandingEllipsoid->constraint("ellipsoid", Expr::vstack(this->alphaExpandingEllipsoid, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, this->xExpandingEllipsoid))), Domain::inQCone());
	}
	else
	{
		this->MExpandingEllipsoid->getConstraint("ellipsoid")->update(Expr::vstack(this->alphaExpandingEllipsoid, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, this->xExpandingEllipsoid))));
	}
	this->MExpandingEllipsoid->solve();
	p_out = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(xExpandingEllipsoid->level()->begin(), n));
	return (*this->alphaExpandingEllipsoid->level())[0];
}

double Circle::closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out)
{
	Eigen::VectorXd v=p_in - this->centroid;
	double vNorm = v.norm();
	if (vNorm < this->radius)
	{
		p_out = p_in;
		double distance = 0.0;
		return distance;
	}
	v = v / vNorm;
	p_out = v * this->radius+this->centroid;
	double distance=(p_out - p_in).norm();
	return distance;
}

double Circle::closestPoint(const Eigen::VectorXd& p_in)
{
	Eigen::VectorXd p_out;
	return this->closestPoint(p_in, p_out);
}

bool Circle::isInside(const Eigen::VectorXd& p, const double &tol)
{
	return (((this->centroid-p)).norm() <= (this->radius+tol));
}

bool Circle::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol)
{
	return ((ai.dot(this->centroid) - bi - this->radius) <=tol);
}

Eigen::VectorXd Circle::getCentroid()
{
	return this->centroid;
}

Eigen::MatrixXd Circle::getBoundingBox()
{
	if (!this->bbComputed)
	{
		this->bb.resize(this->n, 2);
		Eigen::VectorXd q = this->radius * Eigen::VectorXd(this->n);
		this->bb.col(0) = this->centroid - q;
		this->bb.col(1) = this->centroid + q;
		this->bbComputed = true;
	}
	return this->bb;
}

void Circle::allocateClosestPointEllipsoidSolver()
{
	if (!this->solverExpandingEllipsoidAllocated)
	{
		Range* range = Range::getInstance();
		std::pair<double, double> bounds = range->getBounds();
		this->MExpandingEllipsoid = new Model("Circle::closestPointExpandingEllipsoid");
		this->xExpandingEllipsoid = this->MExpandingEllipsoid->variable("x", n, Domain::inRange(bounds.first, bounds.second));
		this->alphaExpandingEllipsoid = this->MExpandingEllipsoid->variable("alpha", 1, Domain::greaterThan(0.0));
		this->MExpandingEllipsoid->constraint("circle", Expr::vstack(this->radius,this->xExpandingEllipsoid), Domain::inQCone());
		this->MExpandingEllipsoid->objective(ObjectiveSense::Minimize, this->alphaExpandingEllipsoid);
		this->solverExpandingEllipsoidAllocated = true;
	}
}
