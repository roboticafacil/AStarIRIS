#include "Sphere.h"
#include "Ellipsoid.h"
#include "ConicSet.h"
#include "Range.h"
#include "EigenNdArray.h"


Sphere::Sphere(const int& n): radius(0.), ConicSet(n) //ConicSet(n,false)
{
}

Sphere::Sphere(const Sphere& sphere): radius(sphere.radius), ConicSet(sphere)
{
}

Sphere::Sphere(const Eigen::VectorXd& center, const double& radius): radius(radius), ConicSet(center.rows()) //ConicSet(center.rows(),false)
{
	this->centroid = center;
}

Sphere::~Sphere()
{
	//auto _M = finally([&]() { M->dispose(); });
	if (solverExpandingEllipsoidAllocated)
	{
		this->MExpandingEllipsoid->dispose();
		solverExpandingEllipsoidAllocated = false;
	}
}

Sphere& Sphere::operator=(const Sphere& other) {
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

double Sphere::closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out)
{
	std::shared_ptr<ndarray<double, 2>> Cinv_ptr = Eigen2NdArray(ellipsoid.getCInv());
	std::shared_ptr<ndarray<double, 1>> centroid_ptr = Eigen2NdArray(ellipsoid.getCentroid());
	this->allocateClosestPointEllipsoidSolver();
	if (!this->MExpandingEllipsoid->hasConstraint("sphere"))
	{
		this->MExpandingEllipsoid->constraint("sphere", Expr::vstack(this->alphaExpandingEllipsoid, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, this->xExpandingEllipsoid))), Domain::inQCone());
	}
	else
	{
		this->MExpandingEllipsoid->getConstraint("sphere")->update(Expr::vstack(this->alphaExpandingEllipsoid, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, this->xExpandingEllipsoid))));
	}
	this->MExpandingEllipsoid->solve();
	p_out = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(xExpandingEllipsoid->level()->begin(), n));
	return (*this->alphaExpandingEllipsoid->level())[0];
}

double Sphere::closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out)
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

double Sphere::closestPoint(const Eigen::VectorXd& p_in)
{
	Eigen::VectorXd p_out;
	return this->closestPoint(p_in, p_out);
}

bool Sphere::isInside(const Eigen::VectorXd& p, const double &tol)
{
	return (((this->centroid-p)).norm() <= (this->radius+tol));
}

bool Sphere::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol)
{
	return ((ai.dot(this->centroid) - bi - this->radius) <=tol);
}

bool Sphere::isInsideSeparatingHyperplanes(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const double& tol)
{
	std::cout << "Sphere: isInsideSeparatingHyperplanes not implemented yet!!" << std::endl;
	//Hay que comprobar esto...
	return ((((A * this->centroid) - b).array()) <= (tol+this->radius)).all();
}

Eigen::VectorXd Sphere::getCentroid()
{
	return this->centroid;
}

Eigen::MatrixXd Sphere::getBoundingBox()
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

void Sphere::allocateClosestPointEllipsoidSolver()
{
	if (!this->solverExpandingEllipsoidAllocated)
	{
		Range* range = Range::getInstance();
		std::pair<double, double> bounds = range->getBounds();
		this->MExpandingEllipsoid = new Model("Sphere::closestPointExpandingEllipsoid");
		this->xExpandingEllipsoid = this->MExpandingEllipsoid->variable("x", n, Domain::inRange(bounds.first, bounds.second));
		this->alphaExpandingEllipsoid = this->MExpandingEllipsoid->variable("alpha", 1, Domain::greaterThan(0.0));
		this->MExpandingEllipsoid->constraint("sphere", Expr::vstack(this->radius,this->xExpandingEllipsoid), Domain::inQCone());
		this->MExpandingEllipsoid->objective(ObjectiveSense::Minimize, this->alphaExpandingEllipsoid);
		this->solverExpandingEllipsoidAllocated = true;
	}
}
