#include<Eigen/Dense>
#include "Polyhedron.h"
//#include "Range.h"
#include "fusion.h"
#include "EigenNdArray.h"
#include "Ellipsoid.h"
#include "Range.h"

using namespace orgQhull;

Polyhedron::Polyhedron(const int& n): ConicSet(n), ellipsoid(n), circle(n) //ConicSet(n,true)
{
}

/*Polyhedron::Polyhedron(const Eigen::MatrixXd& v, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Range& range, const bool allocateSolver) : v(v), A(A), b(b), ConicSet(v.rows()) //ConicSet(v.rows(),true)
{
    if (v.cols() > 0)
    {
        this->v = v;
        if (A.cols() > 0)
        {
            assert(A.rows() == b.rows());
        }
        else
        {
            Eigen::MatrixXd tmpA(A.rows(), A.cols());
            Eigen::VectorXd tmpB(b.rows());
            Polyhedron::vert2con(this->v, tmpA, tmpB);
            this->A=Eigen::MatrixXd(A.rows() + range.A.rows(), A.cols());
            this->b=Eigen::VectorXd(b.rows() + range.b.rows());
            this->A << tmpA, range.A;
            this->b << tmpB, range.b;
        }
    }
    else
    {
        this->v = v;
        if (A.cols() > 0)
        {
            assert(A.rows() == b.rows());
        }
        this->A = Eigen::MatrixXd(A.rows() + range.A.rows(), A.cols());
        this->b = Eigen::VectorXd(b.rows() + range.b.rows());
        this->A << A, range.A;
        this->b << b, range.b;
    }
    //std::cout << this->A << std::endl;
    A_ptr = Eigen2NdArray(this->A);
    b_ptr = Eigen2NdArray(this->b);
    if (allocateSolver)
    {
        //Prepare optimization for closest point computation (this optimization problem is supposed to be called several times, so it's worthy to allocate the variables and constraints in the constructor)
        MclosestPoint = new Model("Polyhedron::closestPoint");
        xclosestPoint = MclosestPoint->variable("x", n, Domain::unbounded());
        l = MclosestPoint->variable("l", 1, Domain::greaterThan(0.0));
        MclosestPoint->constraint("linear", Expr::sub(Expr::mul(A_ptr, xclosestPoint), b_ptr), Domain::lessThan(0.0));
        MclosestPoint->objective(ObjectiveSense::Minimize, l);
        this->solverAllocated = allocateSolver;
    }
}*/

Polyhedron::Polyhedron(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) : A(A), b(b), ConicSet(A.cols()), ellipsoid(A.cols()), circle(A.cols()) //ConicSet(v.rows(),true)
{
    A_ptr = Eigen2NdArray(this->A);
    b_ptr = Eigen2NdArray(this->b);
}
Polyhedron::Polyhedron(const Polyhedron& polyhedron): ConicSet(polyhedron), A(polyhedron.A), b(polyhedron.b), ellipsoid(polyhedron.ellipsoid), circle(polyhedron.circle)
{
    A_ptr = Eigen2NdArray(this->A);
    b_ptr = Eigen2NdArray(this->b);
}
Polyhedron::~Polyhedron()
{
    //auto _MclosestPoint = finally([&]() { MclosestPoint->dispose(); });
}
void Polyhedron::update(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
{
    int oldSize = this->A.rows();
    int newSize = A.rows();
    this->A = A;
    this->b = b;
    A_ptr = Eigen2NdArray(this->A);
    b_ptr = Eigen2NdArray(this->b);
    //First dispose previously allocated models
    if (oldSize != newSize)
    {
        if (this->solverClosestPointEllipsoidAllocated)
        {
            this->solverClosestPointEllipsoidAllocated = false;
            this->MClosestPointEllipsoid->dispose();
        }
        if (this->solverInscribedEllipsoidAllocated)
        {
            this->solverInscribedEllipsoidAllocated = false;
            this->MInscribedEllipsoid->dispose();
        }
        if (this->solverClosestPointAllocated)
        {
            this->solverClosestPointAllocated = false;
            this->MClosestPoint->dispose();
        }
        if (this->solverIsInsideSeparatingHyperplaneAllocated)
        {
            this->solverIsInsideSeparatingHyperplaneAllocated = false;
            this->MIsInsideSeparatingHyperplane->dispose();
        }
    }
    /*if (this->solverClosestPointEllipsoidAllocated)
    {
        if (this->MClosestPointEllipsoid->hasConstraint("linear"))
        {
            this->MClosestPointEllipsoid->getConstraint("linear")->update(Expr::sub(Expr::mul(A_ptr, this->xClosestPointEllipsoid), b_ptr));
        }
        else
            this->MClosestPointEllipsoid->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->xClosestPointEllipsoid), b_ptr), Domain::lessThan(0.0));
    }
    if (this->solverInscribedEllipsoidAllocated)
    {
        if (!this->MInscribedEllipsoid->hasConstraint("Q"))
        {
            this->MInscribedEllipsoid->constraint("Q", Expr::hstack(Expr::sub(b_ptr, Expr::mul(A_ptr, this->xInscribedEllipsoid)), Expr::mul(A_ptr, this->CInscribedEllipsoid)), Domain::inQCone());
        }
        else
        {
            this->MInscribedEllipsoid->getConstraint("Q")->update(Expr::hstack(Expr::sub(b_ptr, Expr::mul(A_ptr, this->xInscribedEllipsoid)), Expr::mul(A_ptr, this->CInscribedEllipsoid)));
        }
    }
    if (this->solverClosestPointAllocated)
    {
        if (this->MClosestPoint->hasConstraint("linear"))
        {
            this->MClosestPoint->getConstraint("linear")->update(Expr::sub(Expr::mul(A_ptr, this->xClosestPoint), b_ptr));
        }
        else
            this->MClosestPoint->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->xClosestPoint), b_ptr), Domain::lessThan(0.0));
    }
    if (this->solverIsInsideSeparatingHyperplaneAllocated)
    {
        if (this->MIsInsideSeparatingHyperplane->hasConstraint("linear"))
        {
            this->MIsInsideSeparatingHyperplane->getConstraint("linear")->update(Expr::sub(Expr::mul(A_ptr, this->xIsInsideSeparatingHyperplane), b_ptr));
        }
        else
            this->MIsInsideSeparatingHyperplane->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->xIsInsideSeparatingHyperplane), b_ptr), Domain::lessThan(0.0));
    }*/
}
bool Polyhedron::isInside(const Eigen::VectorXd& q)
{
    //std::cout << this->A << std::endl;
    //std::cout << this->b << std::endl;
    //std::cout << q << std::endl;
    //bool res = ((this->A * q - this->b).array() <= 0.).all();
    //std::cout << res << std::endl;
    return ((this->A * q - this->b).array()<=0.).all();
}

bool Polyhedron::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi)
{
    this->allocateIsInsideSeparatingHyperplaneSolver();
    std::shared_ptr<ndarray<double, 2>> ai_ptr = Eigen2NdArray(Eigen::MatrixXd(ai.transpose()));
    if (!this->MIsInsideSeparatingHyperplane->hasConstraint("hyperplane"))
        this->MIsInsideSeparatingHyperplane->constraint("hyperplane", Expr::sub(Expr::mul(ai_ptr, this->xIsInsideSeparatingHyperplane), bi), Domain::lessThan(0.0));
    else
        this->MIsInsideSeparatingHyperplane->getConstraint("hyperplane")->update(Expr::sub(Expr::mul(ai_ptr, this->xIsInsideSeparatingHyperplane), bi));
    this->MIsInsideSeparatingHyperplane->solve();
    ProblemStatus status = this->MIsInsideSeparatingHyperplane->getProblemStatus();
    return (status == ProblemStatus::PrimalFeasible) || (status == ProblemStatus::PrimalAndDualFeasible);
}

void Polyhedron::vert2con(const Eigen::MatrixXd& v, Eigen::MatrixXd& A_out, Eigen::VectorXd& b_out)
{
    const int n = v.rows();
    const int num_points = v.cols();
    Eigen::MatrixXd vcpy = v;
    //std::cout << "vcp=" << std::endl;
    //std::cout << vcpy << std::endl;
    std::string comment = ""; // rbox commands, see http://www.qhull.org/html/rbox.htm
    std::string qhull_command = ""; // For qhull commands, see http://www.qhull.org/html/qhull.htm
    try
    {
        Qhull qhull;
        qhull.runQhull(comment.c_str(), n, num_points, vcpy.data(), qhull_command.c_str());
        //std::cout << "qhull # vertices: " << qhull.vertexList().size() << std::endl;
        //std::cout << "qhull.hullDimension(): " << qhull.hullDimension() << std::endl;
        //std::cout << "qhull.volume(): " << qhull.volume() << std::endl;
        //std::cout << "qhull.area(): " << qhull.area() << std::endl;
        const int num_vertex = qhull.vertexList().size();
        QhullVertexList::iterator it = qhull.vertexList().begin();
        //Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> cc(qhull.vertexList().begin()->getVertexT()->point, num_vertex, n);
        Eigen::MatrixXd cc;
        cc.resize(num_vertex, n);
        for (int i = 0; i < num_vertex; i++)
        {
            double* coords = it->point().coordinates();
            for (int j = 0; j < n; j++)
                cc(i, j) = coords[j];
            it++;
        }
        //std::cout << "cc=" << std::endl;
        //std::cout << cc << std::endl;
        Eigen::VectorXd c;
        c.resize(n);
        for (int i = 0; i < n; i++)
        {
            double m = cc.col(i).mean();
            c(i) = m;
            for (int j = 0; j < num_points; j++)
            {
                vcpy(i, j) -= m;
            }
        }
        Eigen::MatrixXd Atmp;
        Atmp.resize(num_vertex, n);
        //std::cout << vcpy << std::endl;
        //std::cout << "Centroid=" << std::endl;
        //std::cout << c << std::endl;
        QhullFacetList facets = qhull.facetList();
        //std::cout << "Facets: " << facets.size() << std::endl;
        int rc = 0;
        for (QhullFacetList::iterator it = facets.begin(); it != facets.end(); ++it)
        {
            if (!(*it).isGood()) continue;
            QhullFacet f = *it;
            QhullVertexSet vSet = f.vertices();
            //std::cout << "Facet: " << std::endl;
            QhullVertexSet::iterator vIt = vSet.begin();
            //QhullVertex vertex = *vIt;
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> F;
            F.resize(vSet.size(), n);
            int j = 0;
            //
            for (vIt = vSet.begin(); vIt != vSet.end(); ++vIt)
            {
                QhullVertex vertex = *vIt;
                QhullPoint p = vertex.point();
                double* coords = p.coordinates();
                //std::cout << coords[0] << "," << coords[1] << std::endl;
                for (int k = 0; k < n; k++)
                    F(j, k) = coords[k];
                j++;
            }
            //std::cout << F << std::endl;
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp(F);
            Eigen::Index r = decomp.rank();

            //std::cout << "Rank: ";
            //std::cout << r << std::endl;
            if (r == F.rows())
            {
                Eigen::VectorXd calc = decomp.solve(Eigen::VectorXd::Ones(F.rows()));
                //std::cout << calc.transpose() << std::endl;
                Atmp.row(rc) = calc.transpose();
                rc++;
            }
            //
            //std::cout << Atmp << std::endl;
            /*for (QhullVertexSet::iterator vIt = vSet.begin(); vIt != vSet.end(); ++vIt)
            {
                QhullVertex vertex = *vIt;
                QhullPoint p = vertex.point();
                double* coords = p.coordinates();
                std::cout << coords[0] << "," << coords[1] << std::endl;
                //vec3 aPoint = vec3(coords[0], coords[1], coords[2]);
                // ...Do what ever you want
            }*/
        }
        Atmp = Atmp(Eigen::seq(0, rc - 1), Eigen::all);
        //std::cout << "Atmp = " << std::endl;
        //std::cout << Atmp << std::endl;
        Eigen::VectorXd btmp = Eigen::VectorXd::Ones(Atmp.rows());
        btmp += Atmp * c;
        //std::cout << "btmp = " << std::endl;
        //std::cout << btmp << std::endl;
        //A_out = Atmp;
        //b_out = btmp;
        removeRepeatedConstraints(Atmp, btmp, A_out, b_out);
        //std::cout << "A_out=" << std::endl;
        //std::cout << A_out << std::endl;
        //std::cout << "b_out=" << std::endl;
        //std::cout << b_out << std::endl;
    }
    catch (QhullError& e)
    {
        std::cerr << e.what() << std::endl;
    }
}

void Polyhedron::removeRepeatedConstraints(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& A_out, Eigen::VectorXd& b_out)
{
    A_out = A;
    b_out = b;
    //By normalizing and ensuring that the b coefficients are positive we can distinguish between repeated constraints
    for (int i = 0; i < A_out.rows(); i++)
    {
        double nn = A_out.row(i).norm();
        A_out.row(i) /= nn;
        b_out(i) /= nn;
        /*if (b_out(i) < 0)
        {
            b_out(i) = -b_out(i);
            A_out.row(i) = -A_out.row(i);
        }*/
    }
    //std::cout << "A_out=" << std::endl;
    //std::cout << A_out << std::endl;
    //std::cout << "b_out=" << std::endl;
    //std::cout << b_out << std::endl;

    Eigen::MatrixXd Ab(A_out.rows(), A_out.cols() + 1);
    Ab << A_out, b_out;


    std::vector<Eigen::VectorXd> vec;
    for (int i = 0; i < Ab.rows(); ++i)
        vec.push_back(Ab.row(i));

    //std::cout << "Ab=" << std::endl;
    //std::cout << Ab << std::endl;

    std::sort(vec.begin(), vec.end(), [](Eigen::VectorXd const& t1, Eigen::VectorXd const& t2) {
        bool t1Smaller = true;
        for (int i = 0; i < t1.rows(); i++)
        {
            if (abs(t1(i) - t2(i)) < 1e-4)
                continue;
            else
            {
                t1Smaller = (t1(i) < t2(i));
                break;
            }
        }
        return t1Smaller;
        //return t1(0) < t2(0); 
        });

    auto it = std::unique(vec.begin(), vec.end(), [](const Eigen::VectorXd& t1, const Eigen::VectorXd& t2) {
        //std::cout << "t1=" << t1.transpose() << " t2=" << t2.transpose() << " return=" << ((t1 - t2).cwiseAbs().maxCoeff() < 1e-1) << std::endl;
        return (t1 - t2).cwiseAbs().maxCoeff() < 1e-1;
        });
    vec.erase(it, vec.end());
    //vec.resize(std::distance(vec.begin(), it));

    //std::cout << "vec.size() = " << vec.size() << std::endl;

    A_out.resize(vec.size(), A_out.cols());
    b_out.resize(vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
        for (int j = 0; j < A_out.cols(); j++)
            A_out(i, j) = vec[i](j);
        b_out(i) = vec[i](A_out.cols());
    }
};

void Polyhedron::con2vert(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& v)
{
    //TO BE IMPLEMENTED!!
}

void Polyhedron::print()
{
    std::cout << "A" << std::endl;
    std::cout << this->A << std::endl;
    std::cout << "b=" << std::endl;
    std::cout << this->b << std::endl;
}

double Polyhedron::closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out)
{
    this->allocateClosestPointSolver();
    //Model::t M = new Model("Polyhedron::closestPoint");
    //auto _M = finally([&]() { M->dispose(); });
    //Variable::t x = MclosestPoint->variable("x", n, Domain::unbounded());
    //Point
    std::shared_ptr<ndarray<double, 1>> p_ptr = Eigen2NdArray(p_in);
    //Variable::t l= MclosestPoint->variable("l", 1, Domain::greaterThan(0.0));
    //Constraint::t Axb= MclosestPoint->constraint(Expr::sub(Expr::mul(A_ptr, x), b_ptr), Domain::lessThan(0.0));
    if (!this->MClosestPoint->hasConstraint("Q"))
        this->MClosestPoint->constraint("Q", Expr::vstack(l, Expr::sub(p_ptr, this->xClosestPoint)), Domain::inQCone());
    else
        this->MClosestPoint->getConstraint("Q")->update(Expr::vstack(this->l, Expr::sub(p_ptr, this->xClosestPoint)));
    //Variable::t x = M->variable("x", n, Domain::unbounded());
    //M->constraint(Expr::vstack(l,Expr::sub(p_ptr,x)), Domain::inQCone());    
    this->MClosestPoint->solve();
    //SolutionStatus sol = Mcpy->getPrimalSolutionStatus();
    //std::cout << sol << std::endl;
    p_out=Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->xClosestPoint->level()->begin(), n));
    Eigen::VectorXd vg =(p_out-p_in);
    double distance = vg.norm();
    //Remove variables (and constraints)
    //l->remove();
    //Axb->remove();
    return distance;
}

double Polyhedron::closestPoint(const Eigen::VectorXd& p_in)
{
    Eigen::VectorXd p_out(p_in.size());
    return this->closestPoint(p_in, p_out);
}

double Polyhedron::closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out)
{
    std::shared_ptr<ndarray<double, 2>> Cinv_ptr = Eigen2NdArray(ellipsoid.getCInv());
    std::shared_ptr<ndarray<double, 1>> centroid_ptr = Eigen2NdArray(ellipsoid.getCentroid());
    this->allocateClosestPointEllipsoidSolver();

    if (!this->MClosestPointEllipsoid->hasConstraint("ellipsoid"))
    {
        this->MClosestPointEllipsoid->constraint("ellipsoid", Expr::vstack(this->alphaClosestPointEllipsoid, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, this->xClosestPointEllipsoid))), Domain::inQCone());
    }
    else
    {
        this->MClosestPointEllipsoid->getConstraint("ellipsoid")->update(Expr::vstack(this->alphaClosestPointEllipsoid, Expr::mul(Cinv_ptr, Expr::sub(centroid_ptr, this->xClosestPointEllipsoid))));
    }
    this->MClosestPointEllipsoid->solve();
    p_out = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(xClosestPointEllipsoid->level()->begin(), n));
    return (*this->alphaClosestPointEllipsoid->level())[0];
}

Ellipsoid Polyhedron::inscribedEllipsoid()
{
    if (!computedInscribedEllipsoid)
    {
        this->computeInscribedEllipsoid();
    }
    return this->ellipsoid;
}
void Polyhedron::computeInscribedEllipsoid()
{
    this->allocateInscribedEllipsoidSolver();
    /*if (!this->MInscribedEllipsoid->hasConstraint("Q"))
    {
        this->MInscribedEllipsoid->constraint("Q", Expr::hstack(Expr::sub(b_ptr, Expr::mul(A_ptr, this->xInscribedEllipsoid)), Expr::mul(A_ptr,this->CInscribedEllipsoid)), Domain::inQCone());
    }
    else
    {
        this->MInscribedEllipsoid->getConstraint("Q")->update(Expr::hstack(Expr::sub(b_ptr, Expr::mul(A_ptr, this->xInscribedEllipsoid)), Expr::mul(A_ptr,this->CInscribedEllipsoid)));
    }*/
    this->MInscribedEllipsoid->solve();
    Eigen::MatrixXd CEigen = Eigen::MatrixXd(Eigen::Map<Eigen::MatrixXd>(this->CInscribedEllipsoid->level()->begin(), n, n));
    Eigen::VectorXd dEigen = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->xInscribedEllipsoid->level()->begin(), n));
    this->ellipsoid= Ellipsoid(CEigen, dEigen);
}

Circle Polyhedron::inscribedCircle()
{
    if (!computedInscribedCircle)
    {
        this->computeInscribedCircle();
    }
    return this->circle;
}

void Polyhedron::computeInscribedCircle()
{
    //std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(polyhedron.A);
    //std::shared_ptr<ndarray<double, 1>> b_ptr = Eigen2NdArray(polyhedron.b);
    //Model::t M = new Model("inscribed_ellipsoid"); auto _M = finally([&]() { M->dispose(); });
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    Model::t M = new Model("Polyhedron::inscribedCircle");
    auto _M = finally([&]() { M->dispose(); });
    Variable::t x = M->variable("x", n, Domain::inRange(bounds.first, bounds.second));
    // Setup variables
    Variable::t c = M->variable("c", n, Domain::unbounded());
    Variable::t r = M->variable("r", 1, Domain::greaterThan(0.0));
    Constraint::t Axb=M->constraint(Expr::sub(Expr::mul(A_ptr, x), b_ptr), Domain::lessThan(0.0));
    M->constraint(Expr::vstack(r, Expr::sub(c, x)), Domain::inQCone());
    M->objective(ObjectiveSense::Maximize, r);
    M->solve();
    Eigen::VectorXd cEigen = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(c->level()->begin(), n));
    double radius=(*r->level())[0];
    this->circle=Circle(cEigen, radius);
}

Eigen::VectorXd Polyhedron::getCentroid()
{
    if (!this->centroidComputed)
    {
        Circle circle = this->inscribedCircle();
        this->centroidComputed = true;
        this->centroid = circle.getCentroid();
    }
    return this->centroid;
}

bool Polyhedron::hasVertices()
{
    return false;
}

Polyhedron Polyhedron::intersection(const Polyhedron& polyhedron1, const Polyhedron& polyhedron2)
{
    assert(polyhedron1.A.cols() == polyhedron2.A.cols());
    Eigen::MatrixXd bigA(polyhedron1.A.rows()+polyhedron2.A.rows(),polyhedron1.A.cols());
    Eigen::VectorXd bigB(bigA.rows());
    bigA << polyhedron1.A, polyhedron2.A;
    bigB << polyhedron1.b, polyhedron2.b;
    return Polyhedron(bigA,bigB);
}

void Polyhedron::allocateClosestPointSolver()
{
    if (!this->solverClosestPointAllocated)
    {
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        this->MClosestPoint = new Model("Polyhedron::closestPoint");
        this->xClosestPoint = this->MClosestPoint->variable("x", n, Domain::inRange(bounds.first, bounds.second));
        this->l = this->MClosestPoint->variable("l", 1, Domain::greaterThan(0.0));
        if (this->A.rows() > 0)
        {
            this->MClosestPoint->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->xClosestPoint), b_ptr), Domain::lessThan(0.0));
        }
        this->MClosestPoint->objective(ObjectiveSense::Minimize, this->l);
        this->solverClosestPointAllocated = true;
    }
}

void Polyhedron::allocateClosestPointEllipsoidSolver()
{
    if (!this->solverClosestPointEllipsoidAllocated)
    {
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        this->MClosestPointEllipsoid = new Model("Polyhedron::closestPointEllipsoid");
        this->xClosestPointEllipsoid = this->MClosestPointEllipsoid->variable("x", n, Domain::inRange(bounds.first, bounds.second));
        this->alphaClosestPointEllipsoid = this->MClosestPointEllipsoid->variable("alpha", 1, Domain::greaterThan(0.0));
        if (this->A.rows() > 0)
        {
            this->MClosestPointEllipsoid->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->xClosestPointEllipsoid), b_ptr), Domain::lessThan(0.0));
        }
        this->MClosestPointEllipsoid->objective(ObjectiveSense::Minimize, this->alphaClosestPointEllipsoid);
        this->solverClosestPointEllipsoidAllocated = true;
    }
}

void Polyhedron::allocateIsInsideSeparatingHyperplaneSolver()
{
    if (!this->solverIsInsideSeparatingHyperplaneAllocated)
    {
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        this->MIsInsideSeparatingHyperplane = new Model("Polyhedron::separatingHyperplane");
        this->xIsInsideSeparatingHyperplane = this->MIsInsideSeparatingHyperplane->variable("x", n, Domain::inRange(bounds.first, bounds.second));
        if (this->A.rows() > 0)
        {
            this->MIsInsideSeparatingHyperplane->constraint("linear", Expr::sub(Expr::mul(A_ptr, this->xIsInsideSeparatingHyperplane), b_ptr), Domain::lessThan(0.0));
        }
        this->MIsInsideSeparatingHyperplane->objective(ObjectiveSense::Minimize, Expr::sum(this->xIsInsideSeparatingHyperplane));
        this->solverIsInsideSeparatingHyperplaneAllocated = true;
    }
}

void Polyhedron::allocateInscribedEllipsoidSolver()
{
    if (!this->solverInscribedEllipsoidAllocated)
    {
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        std::pair<std::shared_ptr<ndarray<double, 2>>, std::shared_ptr<ndarray<double, 1>>> rangeConstraints = range->getConstraints();
        this->MInscribedEllipsoid = new Model("Polyhedron::inscribedEllipsoid");
        this->xInscribedEllipsoid = this->MInscribedEllipsoid->variable("x", n, Domain::inRange(bounds.first, bounds.second));
        // Setup variables
        this->tInscribedEllipsoid = this->MInscribedEllipsoid->variable("t", 1, Domain::greaterThan(0.0));
        //Variable::t C = det_rootn(M, t, n);
        //Det rootn
        this->YInscribedEllipsoid = this->MInscribedEllipsoid->variable(Domain::inPSDCone(2 * n));
        this->CInscribedEllipsoid = this->YInscribedEllipsoid->slice(new_array_ptr<int, 1>({ 0, 0 }), new_array_ptr<int, 1>({ n, n }));
        this->ZInscribedEllipsoid = this->YInscribedEllipsoid->slice(new_array_ptr<int, 1>({ 0, n }), new_array_ptr<int, 1>({ n, 2 * n }));
        this->DZInscribedEllipsoid = this->YInscribedEllipsoid->slice(new_array_ptr<int, 1>({ n, n }), new_array_ptr<int, 1>({ 2 * n, 2 * n }));
        // Z is lower-triangular
        std::shared_ptr<ndarray<int, 2>> low_tri(new ndarray<int, 2>(shape_t<2>(n * (n - 1) / 2, 2)));
        int k = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                (*low_tri)(k, 0) = i, (*low_tri)(k, 1) = j, ++k;
        this->MInscribedEllipsoid->constraint(this->ZInscribedEllipsoid->pick(low_tri), Domain::equalsTo(0.0));
        // DZ = Diag(Z)
        this->MInscribedEllipsoid->constraint(Expr::sub(this->DZInscribedEllipsoid, Expr::mulElm(this->ZInscribedEllipsoid, Matrix::eye(n))), Domain::equalsTo(0.0));
        // (Z11*Z22*...*Znn) >= t^n
        this->MInscribedEllipsoid->constraint(Expr::vstack(this->DZInscribedEllipsoid->diag(),this->tInscribedEllipsoid), Domain::inPGeoMeanCone());
        // quadratic cones
        if (this->A.rows() > 0)
        {
            this->MInscribedEllipsoid->constraint("Q", Expr::hstack(Expr::sub(b_ptr, Expr::mul(A_ptr, this->xInscribedEllipsoid)), Expr::mul(A_ptr, this->CInscribedEllipsoid)), Domain::inQCone());
        }
        this->MInscribedEllipsoid->constraint(Expr::hstack(Expr::sub(rangeConstraints.second, Expr::mul(rangeConstraints.first, this->xInscribedEllipsoid)), Expr::mul(rangeConstraints.first, this->CInscribedEllipsoid)), Domain::inQCone());
        // Objective: Maximize t
        this->MInscribedEllipsoid->objective(ObjectiveSense::Maximize,this->tInscribedEllipsoid);
        this->solverInscribedEllipsoidAllocated = true;
    }
}

Eigen::MatrixXd Polyhedron::getBoundingBox()
{
    if (!this->bbComputed)
    {
        this->bb.resize(this->n, 2);
        Model::t bbModel = new Model("Polyhedron::boundinbBox"); auto _M = finally([&]() { bbModel->dispose(); });
        Range* range = Range::getInstance();
        std::pair<double, double> bounds = range->getBounds();
        Variable::t x = bbModel->variable("x", this->n, Domain::inRange(bounds.first, bounds.second));
        bbModel->constraint(Expr::sub(Expr::mul(A_ptr, x), b_ptr), Domain::lessThan(0.0));
        for (int i = 0; i < this->n; i++)
        {
            bbModel->objective(ObjectiveSense::Minimize, x->index(i));
            bbModel->solve();
            bb(i, 0) = bbModel->primalObjValue();
            bbModel->objective(ObjectiveSense::Maximize, x->index(i));
            bbModel->solve();
            bb(i, 1) = bbModel->primalObjValue();
        }
        this->bbComputed = true;
    }
    return this->bb;
}
