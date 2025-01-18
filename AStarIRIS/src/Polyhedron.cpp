#include<Eigen/Dense>
#include "Polyhedron.h"
//#include "Range.h"
#include "fusion.h"
#include "EigenNdArray.h"
#include "Ellipsoid.h"
#include "Range.h"
#include <numbers>

using namespace orgQhull;

std::ostream& operator<<(std::ostream& out, Polyhedron const& p)
{
    {
        out << "A=[";
        for (int i = 0; i < p.A.rows(); i++)
        {
            for (int j = 0; j < p.A.cols(); j++)
            {
                out << p.A(i, j) << " ";
            }
            if (i < (p.A.rows() - 1))
                out << ";";
        }
        out << "];" << std::endl;
        out << "b=[";
        for (int i = 0; i < p.b.rows(); i++)
            out << p.b(i) << " ";
        out << "];" << std::endl;
        return out;
    }
};

Polyhedron::Polyhedron(const int& n): ConicSet(n), ellipsoid(n), sphere(n) //ConicSet(n,true)
{
}

Polyhedron::Polyhedron(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) : A(A), b(b), ConicSet(A.cols()), ellipsoid(A.cols()), sphere(A.cols()) //ConicSet(v.rows(),true)
{
    A_ptr = Eigen2NdArray(this->A);
    b_ptr = Eigen2NdArray(this->b);
}
Polyhedron::Polyhedron(const Polyhedron& polyhedron): ConicSet(polyhedron), A(polyhedron.A), b(polyhedron.b), ellipsoid(polyhedron.ellipsoid), sphere(polyhedron.sphere)
{
    A_ptr = Eigen2NdArray(this->A);
    b_ptr = Eigen2NdArray(this->b);
}

Polyhedron& Polyhedron::operator=(const Polyhedron& other) {
    this->A = other.A;
    this->b = other.b;
    this->A_ptr = other.A_ptr;
    this->b_ptr = other.b_ptr;
    this->bb = other.bb;
    this->bbComputed = other.bbComputed;
    this->centroid = other.centroid;
    this->centroidComputed = other.centroidComputed;
    this->CInscribedEllipsoid = other.CInscribedEllipsoid;
    this->computedInscribedEllipsoid = other.computedInscribedEllipsoid;
    this->computedInscribedSphere = other.computedInscribedSphere;
    this->DZInscribedEllipsoid = other.DZInscribedEllipsoid;
    this->ellipsoid = other.ellipsoid;
    this->l = other.l;
    this->MClosestPoint = other.MClosestPoint;
    this->MClosestPointEllipsoid = other.MClosestPointEllipsoid;
    this->MIsInsideSeparatingHyperplane = other.MIsInsideSeparatingHyperplane;
    this->MInscribedEllipsoid = other.MInscribedEllipsoid;
    this->numIsInsideSeparatingHyperplanes = other.numIsInsideSeparatingHyperplanes;
    this->solverClosestPointAllocated = other.solverClosestPointAllocated;
    this->solverClosestPointEllipsoidAllocated = other.solverClosestPointEllipsoidAllocated;
    this->solverInscribedEllipsoidAllocated = other.solverInscribedEllipsoidAllocated;
    this->solverIsInsideSeparatingHyperplaneAllocated = other.solverIsInsideSeparatingHyperplaneAllocated;
    this->sphere = other.sphere;
    this->tInscribedEllipsoid = other.tInscribedEllipsoid;
    this->wClosestPointEllipsoid = other.wClosestPointEllipsoid;
    this->xClosestPoint = other.xClosestPoint;
    this->xClosestPointEllipsoid = other.xClosestPointEllipsoid;
    this->xInscribedEllipsoid = other.xInscribedEllipsoid;
    this->xIsInsideSeparatingHyperplane = other.xIsInsideSeparatingHyperplane;
    this->YInscribedEllipsoid = other.YInscribedEllipsoid;
    this->ZInscribedEllipsoid = other.ZInscribedEllipsoid;
    return *this;


}
Polyhedron::~Polyhedron()
{
    //auto _MclosestPoint = finally([&]() { MclosestPoint->dispose(); });
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
}
bool Polyhedron::isInside(const Eigen::VectorXd& q, const double &tol)
{
    return ((this->A * q - this->b).array()<=tol).all();
}


bool Polyhedron::isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double& tol)
{
    this->allocateIsInsideSeparatingHyperplaneSolver();
    //if ((ai.rows() != 2))
    //    std::cout << "Something is wrong" << std::endl;
    std::shared_ptr<ndarray<double, 2>> ai_ptr = Eigen2NdArray(Eigen::MatrixXd(ai.transpose()));
    //std::shared_ptr<ndarray<double, 2>> ai_ptr = Eigen2NdArray(Eigen::MatrixXd(ai));
    if (!this->MIsInsideSeparatingHyperplane->hasConstraint("hyperplane"))
        this->MIsInsideSeparatingHyperplane->constraint("hyperplane", Expr::sub(Expr::mul(ai_ptr, this->xIsInsideSeparatingHyperplane), bi), Domain::lessThan(tol));
    else
        this->MIsInsideSeparatingHyperplane->getConstraint("hyperplane")->update(Expr::sub(Expr::mul(ai_ptr, this->xIsInsideSeparatingHyperplane), bi));
    this->numIsInsideSeparatingHyperplanes = 1;
    this->MIsInsideSeparatingHyperplane->solve();
    ProblemStatus status = this->MIsInsideSeparatingHyperplane->getProblemStatus();
    return (status == ProblemStatus::PrimalFeasible) || (status == ProblemStatus::PrimalAndDualFeasible);
}

bool Polyhedron::isInsideSeparatingHyperplanes(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const double& tol)
{
    if (((A.rows() - 1) != this->numIsInsideSeparatingHyperplanes)||(A.rows()==0))
    {
        if (this->solverIsInsideSeparatingHyperplaneAllocated)
            this->MIsInsideSeparatingHyperplane->dispose();
        this->solverIsInsideSeparatingHyperplaneAllocated = false;
        if (A.rows() == 0)
        {
            this->numIsInsideSeparatingHyperplanes = 0;
            return true;
        }
    }
    this->allocateIsInsideSeparatingHyperplaneSolver();
    std::shared_ptr<ndarray<double, 2>> ai_ptr = Eigen2NdArray(Eigen::MatrixXd(A.row(A.rows() - 1)));
    double bi = b(b.rows() - 1);
    this->MIsInsideSeparatingHyperplane->constraint(Expr::sub(Expr::mul(ai_ptr, this->xIsInsideSeparatingHyperplane), bi), Domain::lessThan(tol));
    this->numIsInsideSeparatingHyperplanes = A.rows();
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
        A_out.resize(Atmp.rows(), Atmp.cols());
        b_out.resize(btmp.rows());
        for (int i = 0; i < Atmp.rows(); i++)
        {
            double nA = Atmp.row(i).norm();
            A_out.row(i) = Atmp.row(i)/nA;
            b_out(i) = btmp(i)/nA;
        }
        //removeRepeatedConstraints(Atmp, btmp, A_out, b_out);
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

// Function to compare two rows with tolerance due to floating-point precision issues
bool Polyhedron::areRowsEqual(const Eigen::RowVectorXd& row1, const Eigen::RowVectorXd& row2, const double& tolerance) {
    return (row1 - row2).norm() < tolerance; // Compare the norm of the difference between rows
}

// Function to remove duplicate rows from an Eigen::MatrixXd
Eigen::MatrixXd Polyhedron::removeDuplicateRows(const Eigen::MatrixXd& matrix, const double& tolerance) {
    std::vector<Eigen::RowVectorXd> uniqueRows;

    for (int i = 0; i < matrix.rows(); ++i) {
        Eigen::RowVectorXd currentRow = matrix.row(i);
        bool isDuplicate = false;

        // Check if the current row already exists in uniqueRows
        for (const auto& uniqueRow : uniqueRows) {
            if (Polyhedron::areRowsEqual(currentRow, uniqueRow, tolerance)) {
                isDuplicate = true;
                break;
            }
        }

        // If not a duplicate, add it to uniqueRows
        if (!isDuplicate) {
            uniqueRows.push_back(currentRow);
        }
    }

    // Build the new matrix from uniqueRows
    Eigen::MatrixXd result(uniqueRows.size(), matrix.cols());
    for (size_t i = 0; i < uniqueRows.size(); ++i) {
        result.row(i) = uniqueRows[i];
    }
    return result;
}

std::vector<int> Polyhedron::removeConstraints()
{
    Eigen::MatrixXd Anew;
    Eigen::VectorXd bnew;
    Eigen::MatrixXd AB(A.rows(),A.cols()+1);
    AB << this->A, this->b;
    //this->print();
    //std::cout << AB << std::endl;
    Eigen::MatrixXd ABnew = Polyhedron::removeDuplicateRows(AB);
    //std::cout << ABnew << std::endl;
    //std::cout << ABnew.block(0, 0, ABnew.rows(), AB.cols() - 1) << std::endl;
    //std::cout << ABnew.block(0, AB.cols() - 1, ABnew.rows(), 1) << std::endl;
    //std::vector<int> removed = Polyhedron::removeConstraints(this->A, this->b, Anew, bnew);
    std::vector<int> removed = Polyhedron::removeConstraints(ABnew.block(0,0,ABnew.rows(),AB.cols()-1),ABnew.block(0,AB.cols()-1,ABnew.rows(),1), Anew, bnew);
    if (A.rows() != Anew.rows())
    {
        this->update(Anew, bnew);
    }
    return removed;
}

std::vector<int> Polyhedron::removeConstraints(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& Aout, Eigen::VectorXd& bout)
{
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    Model::t M = new Model("Polyhedron::removeConstraints");
    auto _M = finally([&]() { M->dispose(); });
    Variable::t x = M->variable("x", 2, Domain::inRange(bounds.first, bounds.second));
    Constraint::t con;
    Aout = A;
    bout = b;
    std::shared_ptr<ndarray<double, 2>> A_ptr;
    std::shared_ptr<ndarray<double, 1>> b_ptr;
    A_ptr = Eigen2NdArray(Aout);
    b_ptr = Eigen2NdArray(bout);
    con = M->constraint("linear", Expr::sub(Expr::mul(A_ptr, x), b_ptr), Domain::lessThan(0.0));
    M->objective(ObjectiveSense::Minimize, Expr::sum(x));

    std::vector<int> keep;
    std::vector<int> removed;
    for (int i = 0; i < Aout.rows(); i++)
    {
        Aout.row(i) = -Aout.row(i);
        bout(i) = -bout(i);
        A_ptr = Eigen2NdArray(Aout);
        b_ptr = Eigen2NdArray(bout);
        con->update(Expr::sub(Expr::mul(A_ptr, x), b_ptr));
        //M->writeTask("problem.lp");
        M->solve();
        ProblemStatus status = M->getProblemStatus();
        bool feasible = (status == ProblemStatus::PrimalFeasible) || (status == ProblemStatus::PrimalAndDualFeasible);
        if (feasible)
            keep.push_back(i);
        else
            removed.push_back(i);
        //Undo changes
        Aout.row(i) = -Aout.row(i);
        bout(i) = -bout(i);
    }

    Aout.resize(keep.size(), A.cols());
    bout.resize(keep.size());
    for (int i = 0; i < keep.size(); i++)
    {
        Aout.row(i) = A.row(keep[i]);
        bout(i) = b(keep[i]);
    }
    return removed;
};

double Polyhedron::nchoosek(const int n, const int k)
{
    if (k > n / 2)
    {
        return nchoosek(n, n - k);
    }
    else if (k == 1)
    {
        return n;
    }
    else
    {
        double c = 1;
        for (int i = 1; i <= k; i++)
        {
            c *= (((double)n - k + i) / ((double)i));
        }
        return std::round(c);
    }
}

void Polyhedron::nchoosek(const Eigen::VectorXi& V,const int k,Eigen::MatrixXi& U)
{
    if (V.size() == 0)
    {
        U.resize(0, k);
        return;
    }
    assert((V.cols() == 1 || V.rows() == 1) && "V must be a vector");
    U.resize(nchoosek(V.size(), k), k);
    int running_i = 0;
    int running_j = 0;
    Eigen::VectorXi running(k);
    int N = V.size();
    const std::function<void(int, int)> doCombs =
        [&running, &N, &doCombs, &running_i, &running_j, &U, &V](int offset, int k)
    {
        if (k == 0)
        {
            U.row(running_i) = running;
            running_i++;
            return;
        }
        for (int i = offset; i <= N - k; ++i)
        {
            running(running_j) = V(i);
            running_j++;
            doCombs(i + 1, k - 1);
            running_j--;
        }
    };
    doCombs(0, k);
    U = Eigen::MatrixXi(U.array() - 1);
}

void Polyhedron::unique(Eigen::MatrixXd& A)
{
    std::vector<Eigen::VectorXd> vec;
    for (int64_t i = 0; i < A.rows(); ++i)
        vec.push_back(A.row(i));

    std::sort(vec.begin(), vec.end(), [](Eigen::VectorXd const& t1, Eigen::VectorXd const& t2) { return t1(0) < t2(0); });

    auto it = std::unique(vec.begin(), vec.end());
    vec.resize(std::distance(vec.begin(), it));

    A.resize(vec.size(), A.cols());
    for (int64_t i = 0; i < vec.size(); ++i)
        A.row(i) = vec[i];
}

void Polyhedron::con2vert(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& v)
{
    Eigen::MatrixXi Combs;
    Eigen::VectorXi vv(Eigen::VectorXi::LinSpaced(A.rows(), 1, A.rows()));
    Polyhedron::nchoosek(vv, A.cols(),Combs);
    v.resize(A.cols(),0);
    int k = 0;
    for (int i = 0; i < Combs.rows(); i++)
    {
        Eigen::ArrayXi ind(Combs.cols());
        for (int j = 0; j < Combs.cols(); j++)
        {
            ind(j) = Combs(i, j);
        }
        Eigen::MatrixXd Aeq = A(ind,Eigen::placeholders::all);
        Eigen::VectorXd beq = b(ind);
        Eigen::VectorXd x = Aeq.colPivHouseholderQr().solve(beq);
        //std::cout << "i=" << i << std::endl;
        //
        //std::cout << A * x - b << std::endl;
        if (!((A*x-b).array()<= -1e-3).all())
        {
            //std::cout << "x=" << x << std::endl;
            v.conservativeResize(v.rows(),v.cols()+1);
            v.col(k) = x;
            k++;
        }
    }
    //Polyhedron::unique(v);
}

void Polyhedron::print()
{
    std::cout << "A" << std::endl;
    std::cout << this->A << std::endl;
    std::cout << "b=" << std::endl;
    std::cout << this->b << std::endl;
}

std::ostream& Polyhedron::print(std::ostream& out, const std::string &AName, const std::string& bName, const bool &addRangeLimits)
{
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    Eigen::MatrixXd ARange;
    Eigen::VectorXd bRange;
    range->getConstraints(ARange,bRange);
    out << AName << "=[";
    for (int i = 0; i < this->A.rows(); i++)
    {
        for (int j = 0; j < this->A.cols(); j++)
        {
            if (j < (this->A.cols() - 1))
                out << this->A(i, j) << " ";
            else
                out << this->A(i, j);
        }
        if (!addRangeLimits)
        {
            if (i < (this->A.rows() - 1))
                out << ";";
        }
        else
            out << ";";
    }
    if (addRangeLimits)
    {
        for (int i = 0; i < ARange.rows(); i++)
        {
            for (int j = 0; j < ARange.cols(); j++)
            {
                if (j < (ARange.cols() - 1))
                    out << ARange(i, j) << " ";
                else
                    out << ARange(i, j);
            }
            if (i < (ARange.rows() - 1))
                out << ";";
        }
    }
    out << "];" << std::endl;
    out << bName << "=[";
    for (int i = 0; i < this->b.rows(); i++)
    {
        if (!addRangeLimits)
        {
            if (i < (this->b.rows() - 1))
                out << this->b(i) << ";";
            else
                out << this->b(i);
        }
        else
            out << this->b(i) << ";";
    }
    if (addRangeLimits)
    {
        for (int i = 0; i < bRange.rows(); i++)
        {
            if (i < (bRange.rows() - 1))
                out << bRange(i) << ";";
            else
                out << bRange(i);
        }
    }
    out << "];" << std::endl;
    return out;
};

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


double Polyhedron::inscribedEllipsoidVolume(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
{
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    std::pair<std::shared_ptr<ndarray<double, 2>>, std::shared_ptr<ndarray<double, 1>>> rangeConstraints = range->getConstraints();
    int n = A.cols();
    std::shared_ptr<ndarray<double, 2>> A_ptr= Eigen2NdArray(A);
    std::shared_ptr<ndarray<double, 1>> b_ptr= Eigen2NdArray(b);
    Model::t m = new Model("Polyhedron::inscribedEllipsoidVolume");
    auto _M = finally([&]() { m->dispose(); });
    Variable::t x = m->variable("x", n, Domain::inRange(bounds.first, bounds.second));
    // Setup variables
    Variable::t t = m->variable("t", 1, Domain::greaterThan(0.0));
    Variable::t Y = m->variable(Domain::inPSDCone(2 * n));
    Variable::t C = Y->slice(new_array_ptr<int, 1>({ 0, 0 }), new_array_ptr<int, 1>({ n, n }));
    Variable::t Z = Y->slice(new_array_ptr<int, 1>({ 0, n }), new_array_ptr<int, 1>({ n, 2 * n }));
    Variable::t DZ = Y->slice(new_array_ptr<int, 1>({ n, n }), new_array_ptr<int, 1>({ 2 * n, 2 * n }));
    // Z is lower-triangular
    std::shared_ptr<ndarray<int, 2>> low_tri(new ndarray<int, 2>(shape_t<2>(n * (n - 1) / 2, 2)));
    int k = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            (*low_tri)(k, 0) = i, (*low_tri)(k, 1) = j, ++k;
    m->constraint(Z->pick(low_tri), Domain::equalsTo(0.0));
    // DZ = Diag(Z)
    m->constraint(Expr::sub(DZ, Expr::mulElm(Z, Matrix::eye(n))), Domain::equalsTo(0.0));
    m->constraint(Expr::vstack(DZ->diag(), t), Domain::inPGeoMeanCone());
    // quadratic cones
    m->constraint(Expr::hstack(Expr::sub(b_ptr, Expr::mul(A_ptr, x)), Expr::mul(A_ptr,C)), Domain::inQCone());
    m->constraint(Expr::hstack(Expr::sub(rangeConstraints.second, Expr::mul(rangeConstraints.first, x)), Expr::mul(rangeConstraints.first, C)), Domain::inQCone());
    // Objective: Maximize t
    m->objective(ObjectiveSense::Maximize, t);
    m->solve();
    ProblemStatus status = m->getProblemStatus();
    if ((status == ProblemStatus::PrimalFeasible) || (status == ProblemStatus::PrimalAndDualFeasible))
        return (*t->level())[0];
    else
        return 0.0;
}
void Polyhedron::computeInscribedEllipsoid()
{
    this->allocateInscribedEllipsoidSolver();
    if (!this->MInscribedEllipsoid->hasConstraint("Q"))
    {
        this->MInscribedEllipsoid->constraint("Q", Expr::hstack(Expr::sub(b_ptr, Expr::mul(A_ptr, this->xInscribedEllipsoid)), Expr::mul(A_ptr, this->CInscribedEllipsoid)), Domain::inQCone());
    }
    else
    {
        this->MInscribedEllipsoid->getConstraint("Q")->update(Expr::hstack(Expr::sub(b_ptr, Expr::mul(A_ptr, this->xInscribedEllipsoid)), Expr::mul(A_ptr, this->CInscribedEllipsoid)));
    }
    this->MInscribedEllipsoid->solve();
    Eigen::MatrixXd CEigen = Eigen::MatrixXd(Eigen::Map<Eigen::MatrixXd>(this->CInscribedEllipsoid->level()->begin(), n, n));
    Eigen::VectorXd dEigen = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->xInscribedEllipsoid->level()->begin(), n));
    this->ellipsoid= Ellipsoid(CEigen, dEigen);
}

Sphere Polyhedron::inscribedSphere()
{
    if (!computedInscribedSphere)
    {
        this->computeInscribedSphere();
    }
    return this->sphere;
}

void Polyhedron::computeInscribedSphere()
{
    //std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(polyhedron.A);
    //std::shared_ptr<ndarray<double, 1>> b_ptr = Eigen2NdArray(polyhedron.b);
    //Model::t M = new Model("inscribed_ellipsoid"); auto _M = finally([&]() { M->dispose(); });
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    Model::t M = new Model("Polyhedron::inscribedSphere");
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
    this->sphere= Sphere(cEigen, radius);
}

Eigen::VectorXd Polyhedron::getCentroid()
{
    if (!this->centroidComputed)
    {
        Sphere sphere = this->inscribedSphere();
        this->centroidComputed = true;
        this->centroid = sphere.getCentroid();
    }
    return this->centroid;
}


bool Polyhedron::intersect(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
{
    Eigen::MatrixXd bigA(this->A.rows() + A.rows(), this->A.cols());
    Eigen::VectorXd bigB(bigA.rows());
    bigA << this->A, A;
    bigB << this->b, b;
    std::shared_ptr<ndarray<double, 2>> ABig_ptr = Eigen2NdArray(bigA);
    std::shared_ptr<ndarray<double, 1>> bBig_ptr = Eigen2NdArray(bigB);
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    Model::t MInstersect = new Model("Polyhedron::intersect");
    Variable::t xIntersect = MInstersect->variable("x", this->A.cols(), Domain::inRange(bounds.first, bounds.second));
    MInstersect->constraint("linear", Expr::sub(Expr::mul(ABig_ptr, xIntersect), bBig_ptr), Domain::lessThan(0.0));
    MInstersect->objective(ObjectiveSense::Minimize, Expr::sum(xIntersect));
    MInstersect->solve();
    ProblemStatus status = MInstersect->getProblemStatus();
    return (status == ProblemStatus::PrimalFeasible) || (status == ProblemStatus::PrimalAndDualFeasible);
}

bool Polyhedron::intersect(const Polyhedron& polyhedron1, const Polyhedron& polyhedron2)
{
    assert(polyhedron1.A.cols() == polyhedron2.A.cols());
    Eigen::MatrixXd bigA(polyhedron1.A.rows()+polyhedron2.A.rows(),polyhedron1.A.cols());
    Eigen::VectorXd bigB(bigA.rows());
    bigA << polyhedron1.A, polyhedron2.A;
    bigB << polyhedron1.b, polyhedron2.b;
    std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(bigA);
    std::shared_ptr<ndarray<double, 1>> b_ptr = Eigen2NdArray(bigB);
    Range* range = Range::getInstance();
    std::pair<double, double> bounds = range->getBounds();
    Model::t MInstersect = new Model("Polyhedron::intersect");
    Variable::t xIntersect = MInstersect->variable("x", polyhedron1.A.cols(), Domain::inRange(bounds.first, bounds.second));
    MInstersect->constraint("linear", Expr::sub(Expr::mul(A_ptr, xIntersect), b_ptr), Domain::lessThan(0.0));
    MInstersect->objective(ObjectiveSense::Minimize, Expr::sum(xIntersect));
    MInstersect->solve();
    ProblemStatus status = MInstersect->getProblemStatus();
    return (status == ProblemStatus::PrimalFeasible) || (status == ProblemStatus::PrimalAndDualFeasible);
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
        Model::t bbModel = new Model("Polyhedron::boundinbBox"); 
        auto _M = finally([&]() { bbModel->dispose(); });
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

std::vector<int> Polyhedron::eqConstraints(const Eigen::VectorXd& p, const double& tol)
{
    std::vector<int> activeConstraints;
    for (int i = 0; i < this->A.rows(); i++)
    {
        double x = (this->A.row(i) * p - this->b(i));
        if (abs(x)<=tol)
        {
            activeConstraints.push_back(i);
        }
    }
    return activeConstraints;
}

void Polyhedron::closestConstraint(const Eigen::VectorXd& p, Eigen::VectorXd& ai, double& bi)
{
    
    Eigen::VectorXd::Index idx;
    (this->A * p - this->b).array().maxCoeff(&idx);
    ai = -this->A.row(idx);
    bi = -this->b(idx);
}
