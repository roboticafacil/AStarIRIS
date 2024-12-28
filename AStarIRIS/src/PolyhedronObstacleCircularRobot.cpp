#include<Eigen/Dense>
#include "PolyhedronV.h"
#include "ConicSet.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "EigenNdArray.h"
#include "Range.h"

PolyhedronObstacleCircularRobot::PolyhedronObstacleCircularRobot(const Eigen::MatrixXd& v, const double& r) : PolyhedronV(v)//, ConicSet(v.rows())
{
    this->r = r;
}

PolyhedronObstacleCircularRobot::PolyhedronObstacleCircularRobot(const PolyhedronObstacleCircularRobot& polyhedronCircularRobot): PolyhedronV(polyhedronCircularRobot), r(polyhedronCircularRobot.r)
{
}

double PolyhedronObstacleCircularRobot::closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out)
{
    this->allocateClosestPointSolver();
    
    std::shared_ptr<ndarray<double, 1>> p_ptr = Eigen2NdArray(p_in);
    if (!MClosestPoint->hasConstraint("Q"))
        MClosestPoint->constraint("Q", Expr::vstack(l, Expr::sub(p_ptr, xClosestPoint)), Domain::inQCone());
    else
        MClosestPoint->getConstraint("Q")->update(Expr::vstack(l, Expr::sub(p_ptr, xClosestPoint)));
    MClosestPoint->solve();
    //std::cout << MclosestPoint->getPrimalSolutionStatus() << std::endl;
    p_out = Eigen::Map<Eigen::VectorXd>(xClosestPoint->level()->begin(), n);
    Eigen::VectorXd vg = (p_out - p_in);
    double distance = vg.norm();
    return distance;
}

double PolyhedronObstacleCircularRobot::closestPoint(const Eigen::VectorXd& p_in)
{
    Eigen::VectorXd p_out;
    return this->closestPoint(p_in, p_out);
}

bool PolyhedronObstacleCircularRobot::isInside(const Eigen::VectorXd& p)
{
    //First check if point is contained inside the outer polyhedron
    if (((this->A * p - (this->b.array()+this->r).matrix()).array()<=0.).all())
    {
        //Now, check if it is inside the inner polyhedron (i.e.: assuming r=0)
        if (Polyhedron::isInside(p))
            return true;
        else
        {
            //if (this->hasVertices())
            //{
                for (int i = 0; i < v.cols(); i++)
                {
                    if ((this->v.col(i) - p).norm() <= this->r)
                        return true;
                }
                //If we reach this point is because it might be on the sides (between the inner and outer polyhedron)
                return this->closestPoint(p) < 1e-3;
            //}
            //else
            //    return this->closestPoint(p) < 1e-3;
        }
    }
    else
        return false;

}

bool PolyhedronObstacleCircularRobot::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi)
{
    
    double maxVal = (ai.transpose() * this->v - bi*Eigen::VectorXd::Ones(1, this->v.cols())).maxCoeff();
    //Check first if it the outer polyhedron is inside
    if (maxVal<=this->r)
    {
        //Now, check if it the inner polyhedron is inside
        if (maxVal<=0.)
        {
            return true;
        }
        else
        {
            return Polyhedron::isInsideSeparatingHyperplane(ai,bi);
        }
    }
    else
        return false;
}

Ellipsoid PolyhedronObstacleCircularRobot::inscribedEllipsoid()
{
    //std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(polyhedron.A);
    //std::shared_ptr<ndarray<double, 1>> b_ptr = Eigen2NdArray(polyhedron.b);
    //Model::t M = new Model("inscribed_ellipsoid"); auto _M = finally([&]() { M->dispose(); });
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    Model::t M = new Model("PolyhedronObstacleCircularRobot::inscribedEllipsoid");
    auto _M = finally([&]() { M->dispose(); });
    Parameter::t radius = M->parameter("radius", 1);
    radius->setValue(this->r);
    // Setup variables
    Variable::t t = M->variable("t", 1, Domain::greaterThan(0.0));
    //Variable::t C = det_rootn(M, t, n);
    //Det rootn
    Variable::t Y = M->variable(Domain::inPSDCone(2 * n));
    Variable::t C = Y->slice(new_array_ptr<int, 1>({ 0, 0 }), new_array_ptr<int, 1>({ n, n }));
    Variable::t Z = Y->slice(new_array_ptr<int, 1>({ 0, n }), new_array_ptr<int, 1>({ n, 2 * n }));
    Variable::t DZ = Y->slice(new_array_ptr<int, 1>({ n, n }), new_array_ptr<int, 1>({ 2 * n, 2 * n }));
    // Z is lower-triangular
    std::shared_ptr<ndarray<int, 2>> low_tri(new ndarray<int, 2>(shape_t<2>(n * (n - 1) / 2, 2)));
    int k = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            (*low_tri)(k, 0) = i, (*low_tri)(k, 1) = j, ++k;
    M->constraint(Z->pick(low_tri), Domain::equalsTo(0.0));
    // DZ = Diag(Z)
    M->constraint(Expr::sub(DZ, Expr::mulElm(Z, Matrix::eye(n))), Domain::equalsTo(0.0));
    // (Z11*Z22*...*Znn) >= t^n
    M->constraint(Expr::vstack(DZ->diag(), t), Domain::inPGeoMeanCone());

    //Centroid of the ellipsoid
    Variable::t d = M->variable("d", n, Domain::inRange(bounds.first,bounds.second));
    // quadratic cones
    M->constraint(Expr::hstack(Expr::sub(Expr::add(b_ptr, radius), Expr::mul(A_ptr, d)), Expr::mul(A_ptr, C)), Domain::inQCone());
    // Objective: Maximize t
    M->objective(ObjectiveSense::Maximize, t);
    M->solve();
    Eigen::MatrixXd CEigen = Eigen::MatrixXd(Eigen::Map<Eigen::MatrixXd>(C->level()->begin(), n, n));
    Eigen::VectorXd dEigen = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(d->level()->begin(), n));
    return Ellipsoid(CEigen, dEigen);
}

Circle PolyhedronObstacleCircularRobot::inscribedCircle()
{
    //std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(polyhedron.A);
    //std::shared_ptr<ndarray<double, 1>> b_ptr = Eigen2NdArray(polyhedron.b);
    //Model::t M = new Model("inscribed_ellipsoid"); auto _M = finally([&]() { M->dispose(); });
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    Model::t M = new Model("PolyhedronObstacleCircularRobot::inscribedCircle");
    auto _M = finally([&]() { M->dispose(); });
    Variable::t x = M->variable("x", n, Domain::inRange(bounds.first, bounds.second));
    Variable::t chi1 = M->variable("chi1", v.rows(), Domain::inRange(bounds.first, bounds.second));
    Variable::t chi2 = M->variable("chi2", v.rows(), Domain::inRange(bounds.first, bounds.second));
    Parameter::t radius = M->parameter("radius", 1);
    radius->setValue(this->r);
    //Linear lifted constraint
    M->constraint(Expr::sub(x, Expr::add(chi1, chi2)), Domain::equalsTo(0.0));
    //chi1 belongs to the polyhedron
    M->constraint(Expr::sub(b_ptr, Expr::mul(A_ptr, chi1)), Domain::greaterThan(0.0));
    //norm of chi2 is smaller than radius
    M->constraint(Expr::vstack(radius, chi2), Domain::inQCone());
    // Setup variables
    Variable::t c = M->variable("c", n, Domain::inRange(bounds.first, bounds.second));
    Variable::t r = M->variable("r", 1, Domain::greaterThan(0.0));
    M->constraint(Expr::vstack(r, Expr::sub(c, x)), Domain::inQCone());
    M->objective(ObjectiveSense::Maximize, r);
    M->solve();
    Eigen::VectorXd cEigen = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(c->level()->begin(), n));
    double rr = (*r->level())[0];
    return Circle(cEigen, rr);
}

Eigen::VectorXd PolyhedronObstacleCircularRobot::getCentroid()
{
    if (!this->centroidComputed)
    {
        Circle circle = this->inscribedCircle();
        this->centroidComputed = true;
        this->centroid = circle.getCentroid();
    }
    return this->centroid;
}

double PolyhedronObstacleCircularRobot::closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out)
{
    return Polyhedron::closestPointExpandingEllipsoid(ellipsoid, p_out);
}

void PolyhedronObstacleCircularRobot::allocateClosestPointSolver()
{
    if (!solverClosestPointAllocated)
    { 
        Polyhedron::allocateClosestPointSolver();
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        this->chi1 = this->MClosestPoint->variable("chi1", v.rows(), Domain::inRange(bounds.first,bounds.second));
        this->chi2 = this->MClosestPoint->variable("chi2", v.rows(), Domain::inRange(bounds.first, bounds.second));
        //radius = this->MClosestPoint->parameter("radius", 1);
        //radius->setValue(this->r);
        //Linear lifted constraint
        this->MClosestPoint->constraint(Expr::sub(this->xClosestPoint, Expr::add(this->chi1, this->chi2)), Domain::equalsTo(0.0));
        //chi1 belongs to the polyhedron
        if (this->MClosestPoint->hasConstraint("linear"))
            this->MClosestPoint->getConstraint("linear")->update(Expr::sub(Expr::mul(A_ptr, this->chi1), b_ptr));
        else
            this->MClosestPoint->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->chi1), b_ptr), Domain::lessThan(0.0));
        //norm of chi2 is smaller than radius
        this->MClosestPoint->constraint(Expr::vstack(this->r, this->chi2), Domain::inQCone());
    }
    //MclosestPoint->objective(ObjectiveSense::Minimize, l);
}

void PolyhedronObstacleCircularRobot::allocateClosestPointEllipsoidSolver()
{
    if (!solverClosestPointEllipsoidAllocated)
    {
        Polyhedron::allocateClosestPointEllipsoidSolver();
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        this->chi1ClosestPointEllipsoid = this->MClosestPointEllipsoid->variable("chi1", v.rows(), Domain::inRange(bounds.first, bounds.second));
        this->chi2ClosestPointEllipsoid = this->MClosestPointEllipsoid->variable("chi2", v.rows(), Domain::inRange(bounds.first, bounds.second));
        //radius = this->MClosestPointEllipsoid->parameter("radius", 1);
        //radius->setValue(this->r);
        //Linear lifted constraint
        this->MClosestPointEllipsoid->constraint(Expr::sub(this->xClosestPointEllipsoid, Expr::add(this->chi1ClosestPointEllipsoid, this->chi2ClosestPointEllipsoid)), Domain::equalsTo(0.0));
        //chi1 belongs to the polyhedron
        if (this->MClosestPointEllipsoid->hasConstraint("linear"))
            this->MClosestPointEllipsoid->getConstraint("linear")->update(Expr::sub(Expr::mul(A_ptr, this->chi1ClosestPointEllipsoid), b_ptr));
        else
            this->MClosestPointEllipsoid->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->chi1ClosestPointEllipsoid), b_ptr), Domain::lessThan(0.0));
        //norm of chi2 is smaller than radius
        this->MClosestPointEllipsoid->constraint(Expr::vstack(this->r, this->chi2ClosestPointEllipsoid), Domain::inQCone());
    }
}

void PolyhedronObstacleCircularRobot::allocateIsInsideSeparatingHyperplaneSolver()
{
    if (!solverIsInsideSeparatingHyperplaneAllocated)
    {
        Polyhedron::allocateIsInsideSeparatingHyperplaneSolver();
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        this->chi1IsInsideSeparatingHyperplane = this->MIsInsideSeparatingHyperplane->variable("chi1", v.rows(), Domain::inRange(bounds.first, bounds.second));
        this->chi2IsInsideSeparatingHyperplane = this->MIsInsideSeparatingHyperplane->variable("chi2", v.rows(), Domain::inRange(bounds.first, bounds.second));
        //radius = this->MIsInsideSeparatingHyperplane->parameter("radius", 1);
        //radius->setValue(this->r);
        //Linear lifted constraint
        this->MIsInsideSeparatingHyperplane->constraint(Expr::sub(this->xIsInsideSeparatingHyperplane, Expr::add(this->chi1IsInsideSeparatingHyperplane, this->chi2IsInsideSeparatingHyperplane)), Domain::equalsTo(0.0));
        //chi1 belongs to the polyhedron
        if (this->MIsInsideSeparatingHyperplane->hasConstraint("linear"))
            this->MIsInsideSeparatingHyperplane->getConstraint("linear")->update(Expr::sub(Expr::mul(A_ptr, this->chi1IsInsideSeparatingHyperplane), b_ptr));
        else
            this->MIsInsideSeparatingHyperplane->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->chi1IsInsideSeparatingHyperplane), b_ptr), Domain::lessThan(0.0));
        //norm of chi2 is smaller than radius
        this->MIsInsideSeparatingHyperplane->constraint(Expr::vstack(this->r, this->chi2IsInsideSeparatingHyperplane), Domain::inQCone());
    }
}

Eigen::MatrixXd PolyhedronObstacleCircularRobot::getBoundingBox()
{
    if (!this->bbComputed)
    {
        Polyhedron::getBoundingBox();  //Compute the bb of the contained polyhedron 
        for (int i = 0; i < this->n; i++) //Now adds the radius
        {
            this->bb(i, 0) -= this->r;
            this->bb(i, 1) += this->r;
        }
    }
    return this->bb;
}
