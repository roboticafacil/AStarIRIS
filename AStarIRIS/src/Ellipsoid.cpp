#include<Eigen/Dense>
#include "Ellipsoid.h"
#include "EigenNdArray.h"
#include "Range.h"

Ellipsoid::Ellipsoid(const int& n): solverAllocated(false), ConicSet(n) //ConicSet(n,true)
{
}

Ellipsoid::Ellipsoid(const Eigen::MatrixXd& C, const Eigen::VectorXd& centroid): ConicSet(centroid.rows())// ConicSet(centroid.rows(),true)
{
    assert(C.rows() == centroid.rows());
    assert(C.cols() == C.rows());
    this->C = C;
    this->centroid = centroid;
    this->centroidComputed = true;   
    centroid_ptr = Eigen2NdArray(this->centroid);
}

Ellipsoid::Ellipsoid(const Ellipsoid& ellipsoid): ConicSet(ellipsoid), C(ellipsoid.C), Cinv(ellipsoid.Cinv), inverseComputed(ellipsoid.inverseComputed)
{
}

Ellipsoid::~Ellipsoid()
{
    //auto _M = finally([&]() { M->dispose(); });
}

Ellipsoid& Ellipsoid::operator=(const Ellipsoid& other) {
    this->alpha = other.alpha;
    this->alphaExpandingEllipsoid = other.alphaExpandingEllipsoid;
    this->bb = other.bb;
    this->bbComputed = other.bbComputed;
    this->C = other.C;
    this->centroid = other.centroid;
    this->centroidComputed = other.centroidComputed;
    this->centroid_ptr = other.centroid_ptr;
    this->Cinv = other.Cinv;
    this->Cinv_ptr = other.Cinv_ptr;
    this->inverseComputed = other.inverseComputed;
    this->M = other.M;
    this->MExpandingEllipsoid = other.MExpandingEllipsoid;
    this->MIsInsideSeparatingHyperplane = other.MIsInsideSeparatingHyperplane;
    this->solverAllocated = other.solverAllocated;
    this->solverExpandingEllipsoidAllocated = other.solverExpandingEllipsoidAllocated;
    this->solverIsInsideSeparatingHyperplaneAllocated = other.solverIsInsideSeparatingHyperplaneAllocated;
    this->x = other.x;
    this->xExpandingEllipsoid = other.xExpandingEllipsoid;
    this->xIsInsideSeparatingHyperplane = other.xIsInsideSeparatingHyperplane;
    return *this;
}

void Ellipsoid::print()
{
    std::cout << "C" << std::endl;
    std::cout << C << std::endl;
    std::cout << "centroid" << std::endl;
    std::cout << this->getCentroid() << std::endl;
}

double Ellipsoid::closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out)
{
    this->allocateSolver();
    std::shared_ptr<ndarray<double, 1>> p_ptr = Eigen2NdArray(p_in);
    if (!M->hasConstraint("c1"))
        M->constraint("c1", Expr::vstack(alpha, Expr::sub(p_ptr, x)), Domain::inQCone());
    else
        M->getConstraint("c1")->update(Expr::vstack(alpha, Expr::sub(p_ptr, x)));
    M->solve();
    p_out = Eigen::Map<Eigen::VectorXd>(x->level()->begin(), n);
    Eigen::VectorXd vg = (p_out - p_in);
    double distance = vg.norm();
    return distance;
}
double Ellipsoid::closestPoint(const Eigen::VectorXd& p_in)
{
    Eigen::VectorXd p_out(p_in.size());
    return this->closestPoint(p_in, p_out);
}

double Ellipsoid::closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out)
{
    std::shared_ptr<ndarray<double, 2>> CinvExpanding_ptr = Eigen2NdArray(ellipsoid.getCInv());
    std::shared_ptr<ndarray<double, 1>> centroidExpanding_ptr = Eigen2NdArray(ellipsoid.getCentroid());
    this->allocateClosestPointEllipsoidSolver();
    if (!this->MExpandingEllipsoid->hasConstraint("expandingEllipsoid"))
    {
        this->MExpandingEllipsoid->constraint("expandingEllipsoid", Expr::vstack(this->alphaExpandingEllipsoid, Expr::mul(CinvExpanding_ptr, Expr::sub(centroidExpanding_ptr, this->xExpandingEllipsoid))), Domain::inQCone());
    }
    else
    {
        this->MExpandingEllipsoid->getConstraint("expandingEllipsoid")->update(Expr::vstack(this->alphaExpandingEllipsoid, Expr::mul(CinvExpanding_ptr, Expr::sub(centroidExpanding_ptr, this->xExpandingEllipsoid))));
    }
    this->MExpandingEllipsoid->solve();
    p_out = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(xExpandingEllipsoid->level()->begin(), n));
    return (*this->alphaExpandingEllipsoid->level())[0];
}

bool Ellipsoid::isInside(const Eigen::VectorXd& p)
{
    return ((this->getCInv()*(this->centroid - p)).norm()<=1.);
}

bool Ellipsoid::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi)
{
    this->allocateIsInsideSeparatingHyperplaneSolver();
    std::shared_ptr<ndarray<double, 2>> ai_ptr = Eigen2NdArray(Eigen::MatrixXd(ai.transpose()));
    if (!this->MIsInsideSeparatingHyperplane->hasConstraint("hyperplane"))
        this->MIsInsideSeparatingHyperplane->constraint("hyperplane", Expr::sub(Expr::mul(ai_ptr, this->xIsInsideSeparatingHyperplane), bi), Domain::lessThan(0.0));
    else
        this->MIsInsideSeparatingHyperplane->getConstraint("hyperplane")->update(Expr::sub(Expr::mul(ai_ptr, this->xIsInsideSeparatingHyperplane), bi));
    this->MIsInsideSeparatingHyperplane->solve();
    ProblemStatus status = this->MIsInsideSeparatingHyperplane->getProblemStatus();
    //std::cout << status << std::endl;
    //Eigen::VectorXd p = Eigen::Map<Eigen::VectorXd>(this->xIsInsideSeparatingHyperplane->level()->begin(), n);
    //std::cout << p << std::endl;
    return (status == ProblemStatus::PrimalFeasible) || (status == ProblemStatus::PrimalAndDualFeasible);
}

Eigen::VectorXd Ellipsoid::getCentroid()
{
    return this->centroid;    
}

Eigen::MatrixXd Ellipsoid::getCInv()
{
    if (!this->inverseComputed)
    {
        this->Cinv = C.inverse();
        this->inverseComputed = true;
    }
    return this->Cinv;
}

void Ellipsoid::allocateSolver()
{
    if (!solverAllocated)
    {
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        Cinv_ptr = Eigen2NdArray(this->getCInv());
        M = new Model("Ellipsoid");
        x = M->variable("x", n, Domain::inRange(bounds.first,bounds.second));
        alpha = M->variable("alpha", 1, Domain::greaterThan(0.0));
        //one = M->parameter("one", 1);
        //one->setValue(1.0);
        M->constraint(Expr::vstack(1.0, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, x))), Domain::inQCone());
        M->objective(ObjectiveSense::Minimize, alpha);
        this->solverAllocated = true;
    }
}

Eigen::MatrixXd Ellipsoid::getBoundingBox()
{
    if (!this->bbComputed)
    {
        this->bb.resize(this->n, 2);
        Eigen::MatrixXd Q = this->C * this->C;
        Eigen::VectorXd q = Eigen::VectorXd(Q.diagonal().cwiseSqrt());
        this->bb.col(0) = this->centroid - q;
        this->bb.col(1) = this->centroid + q;
        this->bbComputed = true;
    }
    return this->bb;
}

void Ellipsoid::tangentHyperplane(const Eigen::VectorXd& p_in, Eigen::VectorXd& ai, double& bi)
{
    this->getCInv(); //computes Cinv
    ai=2.*this->Cinv*(this->Cinv*(p_in-this->centroid));
    bi = ai.dot(p_in);
    double nn = ai.norm();
    ai = ai / nn;
    bi = bi / nn;
}

void Ellipsoid::allocateIsInsideSeparatingHyperplaneSolver()
{
    if (!this->solverIsInsideSeparatingHyperplaneAllocated)
    {
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        this->MIsInsideSeparatingHyperplane = new Model("Ellipsoid::separatingHyperplane");
        this->xIsInsideSeparatingHyperplane = this->MIsInsideSeparatingHyperplane->variable("x", n, Domain::inRange(bounds.first, bounds.second));
        //this->oneIsInsideSeparatingHyperplane = this->MIsInsideSeparatingHyperplane->parameter("one", 1);
        //this->oneIsInsideSeparatingHyperplane->setValue(1.0);
        Cinv_ptr = Eigen2NdArray(this->getCInv());
        this->MIsInsideSeparatingHyperplane->constraint(Expr::vstack(1.0, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, xIsInsideSeparatingHyperplane))), Domain::inQCone());
        this->MIsInsideSeparatingHyperplane->objective(ObjectiveSense::Minimize, Expr::sum(this->xIsInsideSeparatingHyperplane));
        this->solverIsInsideSeparatingHyperplaneAllocated = true;
    }
}

void Ellipsoid::allocateClosestPointEllipsoidSolver()
{
    if (!this->solverExpandingEllipsoidAllocated)
    {
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        Cinv_ptr = Eigen2NdArray(this->getCInv());
        this->MExpandingEllipsoid = new Model("Ellipsoid::closestPointExpandingEllipsoid");
        this->xExpandingEllipsoid = this->MExpandingEllipsoid->variable("x", n, Domain::inRange(bounds.first, bounds.second));
        this->alphaExpandingEllipsoid = this->MExpandingEllipsoid->variable("alpha", 1, Domain::greaterThan(0.0));
        this->MExpandingEllipsoid->constraint(Expr::vstack(1.0, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, this->xExpandingEllipsoid))), Domain::inQCone());
        this->MExpandingEllipsoid->objective(ObjectiveSense::Minimize, this->alphaExpandingEllipsoid);
        this->solverExpandingEllipsoidAllocated = true;
    }
}
