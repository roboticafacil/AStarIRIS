//, const Eigen::MatrixXd& A = Eigen::Matrix<double, 0, 0>(), const Eigen::VectorXd& b = Eigen::Vector<double, 0>()

#include<Eigen/Dense>
#include "PolyhedronV.h"
//#include "Range.h"
#include "fusion.h"
#include "EigenNdArray.h"
#include "Ellipsoid.h"
#include "Range.h"

using namespace orgQhull;

PolyhedronV::PolyhedronV(const int& n) : Polyhedron(n)
{
}

PolyhedronV::PolyhedronV(const Eigen::MatrixXd& v) : Polyhedron(Eigen::MatrixXd(0,v.rows()), Eigen::VectorXd(0))
{
    this->v = v;
    Polyhedron::vert2con(this->v, this->A, this->b);
    //std::cout << this->A << std::endl;
    A_ptr = Eigen2NdArray(this->A);
    b_ptr = Eigen2NdArray(this->b);
}
PolyhedronV::PolyhedronV(const PolyhedronV& polyhedron) : Polyhedron(polyhedron), v(polyhedron.v)
{
    A_ptr = Eigen2NdArray(this->A);
    b_ptr = Eigen2NdArray(this->b);
}
PolyhedronV::~PolyhedronV()
{
    //auto _MclosestPoint = finally([&]() { MclosestPoint->dispose(); });
}

void PolyhedronV::print()
{
    std::cout << this->v << std::endl;
}

bool PolyhedronV::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi)
{
    return (((ai.transpose() * this->v - bi * Eigen::VectorXd::Ones(1, this->v.cols()))).array() <= 0.).any();
}

double PolyhedronV::closestVertex(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out)
{
    double best_dist = (double)std::numeric_limits<double>::infinity();
    int idx = -1;
    for (int i = 0; i < v.cols(); i++)
    {
        double d = (v.col(i) - p_in).norm();
        if (d < best_dist)
        {
            idx = i;
            d = best_dist;
        }
    }
    p_out = v.col(idx);
    return best_dist;
}

double PolyhedronV::closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out)
{
    Eigen::MatrixXd vg = ellipsoid.getCInv() * (this->v - ellipsoid.getCentroid().replicate(1, v.cols()));
    std::shared_ptr<ndarray<double, 2>> vg_ptr = Eigen2NdArray(vg);
    this->allocateClosestPointEllipsoidSolver();

    if (!this->MClosestPointEllipsoidQP->hasConstraint("eq"))
    {
        this->MClosestPointEllipsoidQP->constraint("eq", Expr::sub(Expr::mul(vg_ptr, this->wClosestPointEllipsoidQP), this->xClosestPointEllipsoidQP), Domain::equalsTo(0.0));
    }
    else
    {
        this->MClosestPointEllipsoidQP->getConstraint("eq")->update(Expr::sub(Expr::mul(vg_ptr, this->wClosestPointEllipsoidQP), this->xClosestPointEllipsoidQP));
    }
    //this->MClosestPointEllipsoidQP->writeTask("closestPointExpandingEllipsoid.ptf");
    this->MClosestPointEllipsoidQP->solve(); 
    //ProblemStatus status = this->MClosestPointEllipsoidQP->getProblemStatus();
    //std::cout << status << std::endl;
    Eigen::VectorXd x = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(xClosestPointEllipsoidQP->level()->begin(), n));
    p_out = ellipsoid.C * x + ellipsoid.getCentroid();
    return (*this->alphaClosestPointEllipsoidQP->level())[0];
}


bool PolyhedronV::hasVertices()
{
    return true;
}

void PolyhedronV::allocateClosestPointEllipsoidSolver()
{
    if (!this->solverClosestPointEllipsoidQPAllocated)
    {
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        this->MClosestPointEllipsoidQP = new Model("Polyhedron::closestPointEllipsoidQP");
        //this->xClosestPointEllipsoidQP = this->MClosestPointEllipsoidQP->variable("x", n, Domain::inRange(bounds.first, bounds.second));
        this->xClosestPointEllipsoidQP = this->MClosestPointEllipsoidQP->variable("x", n, Domain::unbounded()); //Be careful because this "x" is unbounded due to the ball-ellipse space transformation
        this->wClosestPointEllipsoidQP = this->MClosestPointEllipsoidQP->variable("w", v.cols(), Domain::greaterThan(0.0));
        this->alphaClosestPointEllipsoidQP = this->MClosestPointEllipsoidQP->variable("alpha", 1, Domain::greaterThan(0.0));
        this->MClosestPointEllipsoidQP->constraint(Expr::sum(this->wClosestPointEllipsoidQP), Domain::equalsTo(1.0));
        this->MClosestPointEllipsoidQP->constraint(Expr::vstack(this->alphaClosestPointEllipsoidQP, this->xClosestPointEllipsoidQP), Domain::inQCone());
        this->MClosestPointEllipsoidQP->objective(ObjectiveSense::Minimize, this->alphaClosestPointEllipsoidQP);
        this->solverClosestPointEllipsoidQPAllocated = true;
    }
}

Eigen::MatrixXd PolyhedronV::getBoundingBox()
{
    if (!this->bbComputed)
    {
        this->bb.resize(this->n, 2);
        for (int i = 0; i < this->n; i++)
        {
            Eigen::VectorXd r = this->v.row(i);
            bb(i, 0) = r.minCoeff();
            bb(i, 1) = r.maxCoeff();
        }
        this->bbComputed = true;
    }
    return this->bb;
}
