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

using namespace orgQhull;
using namespace mosek::fusion;
using namespace monty;


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
    virtual void print();
    static void vert2con(const Eigen::MatrixXd& v, Eigen::MatrixXd& A_out, Eigen::VectorXd& b_out);
    static void removeRepeatedConstraints(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& A_out, Eigen::VectorXd& b_out);
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
};
#endif