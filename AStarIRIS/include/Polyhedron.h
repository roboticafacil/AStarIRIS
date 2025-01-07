#ifndef POLYHEDRON
#define POLYHEDRON
#include<Eigen/Dense>
#include "libqhullcpp\Qhull.h"
#include "libqhullcpp\QhullError.h"
#include "libqhullcpp\QhullFacet.h"
#include "libqhullcpp\QhullFacetList.h"
#include "libqhullcpp\QhullVertex.h"
#include "libqhullcpp\QhullSet.h"
#include "libqhullcpp\QhullVertexSet.h"
#include "fusion.h"
#include "ConicSet.h"
#include "Ellipsoid.h"
#include "Circle.h"
#include <numbers>
//#ifdef WITH_MATPLOT
//#include <cmath>
//#include <matplot/matplot.h>
//#endif 

using namespace orgQhull;
using namespace mosek::fusion;
using namespace monty;

//#ifdef WITH_MATPLOT
//using namespace matplot;
//#endif

class Polyhedron;

std::ostream& operator<<(std::ostream& out, Polyhedron const& p);

class Polyhedron: public ConicSet
{
public:
    //Eigen::MatrixXd v;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
public:
    Polyhedron(const int &n);
    //Polyhedron(const Eigen::MatrixXd& v, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Range& range, const bool allocateSolver=false);
    Polyhedron(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
    Polyhedron(const Polyhedron& polyhedron);
    ~Polyhedron();
    void update(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
    
    /* {
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
    };*/
    virtual void print();
    std::ostream& print(std::ostream& out, const std::string &AName, const std::string &bName, const bool& addRangeLimits);
    static void vert2con(const Eigen::MatrixXd& v, Eigen::MatrixXd& A_out, Eigen::VectorXd& b_out);
    std::vector<int> removeRepeatedConstraints();
    static std::vector<int> removeRepeatedConstraints(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& A_out, Eigen::VectorXd& b_out);
    static void con2vert(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& v);
    double closestVertex(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    virtual double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    virtual double closestPoint(const Eigen::VectorXd& p_in);
    virtual double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
    virtual bool isInside(const Eigen::VectorXd& p, const double &tol= 1.e-4);
    virtual bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol= 1.e-4);
    virtual Ellipsoid inscribedEllipsoid();
    void computeInscribedEllipsoid();
    virtual Circle inscribedCircle();
    void computeInscribedCircle();
    virtual Eigen::VectorXd getCentroid();
    bool hasVertices();
    static Polyhedron intersection(const Polyhedron& polyhedron1, const Polyhedron& polyhedron2);
    virtual void allocateClosestPointSolver();
    virtual void allocateClosestPointEllipsoidSolver();
    virtual void allocateIsInsideSeparatingHyperplaneSolver();
    void allocateInscribedEllipsoidSolver();
    virtual Eigen::MatrixXd getBoundingBox();
    virtual void getFilled2DPolyhedron(std::vector<double>& x, std::vector<double>& y);
protected:
    //Polyhedron constraints
    std::shared_ptr<ndarray<double, 2>> A_ptr;
    std::shared_ptr<ndarray<double, 1>> b_ptr;
    Model::t MClosestPoint;
    Variable::t xClosestPoint;
    Variable::t l;
    Ellipsoid ellipsoid;
    bool computedInscribedEllipsoid=false;
    Circle circle;
    bool computedInscribedCircle=false;
    bool solverClosestPointAllocated =false;
    Model::t MClosestPointEllipsoid;
    Variable::t xClosestPointEllipsoid;
    Variable::t wClosestPointEllipsoid;
    Variable::t alphaClosestPointEllipsoid;
    bool solverClosestPointEllipsoidAllocated = false;
    Model::t MIsInsideSeparatingHyperplane;
    Variable::t xIsInsideSeparatingHyperplane;
    bool solverIsInsideSeparatingHyperplaneAllocated = false;
    void allocateClosestPointEllipsoidQPSolver();
    Model::t MInscribedEllipsoid;
    Variable::t xInscribedEllipsoid;
    Variable::t CInscribedEllipsoid;
    Variable::t tInscribedEllipsoid;
    Variable::t YInscribedEllipsoid;
    Variable::t ZInscribedEllipsoid;
    Variable::t DZInscribedEllipsoid;
    bool solverInscribedEllipsoidAllocated = false;
private:
    static double nchoosek(const int n, const int k);
    static void nchoosek(const Eigen::VectorXi& V, const int k, Eigen::MatrixXi& U);
    static void unique(Eigen::MatrixXd& A);
};
#endif