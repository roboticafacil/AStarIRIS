#ifndef POLYHEDRONV
#define POLYHEDRONV
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
#include "Polyhedron.h"
#include "Ellipsoid.h"
#include "Circle.h"

using namespace orgQhull;
using namespace mosek::fusion;
using namespace monty;


class PolyhedronV : public Polyhedron
{
public:
    Eigen::MatrixXd v;
public:
    PolyhedronV(const int& n);
    //Polyhedron(const Eigen::MatrixXd& v, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Range& range, const bool allocateSolver=false);
    PolyhedronV(const Eigen::MatrixXd& v);
    PolyhedronV(const PolyhedronV& polyhedron);
    ~PolyhedronV();
    void print();
    double closestVertex(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    virtual double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
    virtual bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol=1.e-4);
    bool hasVertices();
    virtual void allocateClosestPointEllipsoidSolver();
    virtual Eigen::MatrixXd getBoundingBox();
protected:
    Model::t MClosestPointEllipsoidQP;
    Variable::t xClosestPointEllipsoidQP;
    Variable::t wClosestPointEllipsoidQP;
    Variable::t alphaClosestPointEllipsoidQP;
    bool solverClosestPointEllipsoidQPAllocated = false;
};
#endif