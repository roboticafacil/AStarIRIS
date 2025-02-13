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
#include "Sphere.h"
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

/**
  Polyhedron class is used to define a polyhedron \f$\mathcal{P}(A,b)=\{Ax\leq b \forall x\in\mathcal{C}\}\f$, being \f$Ax\leq b\f$ its hyperplanes (the class does not explicitly compute polyhedron vertices) and \f$\mathcal{C}\f$ the configuration space. This class uses Mosek to allocate some solvers for cached performance, i.e.: computing the distance of a point to the polyhedron. As long as the polyhedron doesn't change the number of constraints the cached computation will be kepts in memory.
*/
class Polyhedron: public ConicSet
{
public:
    /**
    Contraints matrix \f$A\in\mathbb{R}^{m\times n}\f$, being \f$m\f$ the number of constraints and \f$n\f$ the polyhedron dimension. Rows of \f$A\f$ must be normalized.
    */
    Eigen::MatrixXd A;
    /**
    Contraints intercepts \f$b\in\mathbb{R}^{m}\f$, being \f$m\f$ the number of constraints.
    */
    Eigen::VectorXd b;
public:
    /**
    Creates an empty polyhedron with the specified dimension.
    \param n Polyhedron dimension
    */
    Polyhedron(const int &n);
    /**
    Creates a polyhedron \f$\mathcal{P}(A,b)=\{Ax\leq b \forall x\in\mathcal{C}\}\f$ with the specified constraints, being \f$\mathcal{C}\f$ the configuration space.
    \param A Contraints matrix \f$A\in\mathbb{R}^{m\times n}\f$, being \f$m\f$ the number of constraints and \f$n\f$ the polyhedron dimension. Rows of \f$A\f$ are normalized
    \param b Contraints intercepts \f$b\in\mathbb{R}^{m}\f$, being \f$m\f$ the number of constraint
    */
    Polyhedron(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
    /**
    Copy constructor.
    \param polyhedron Another polyhedron
    */
    Polyhedron(const Polyhedron& polyhedron);
    /**
    Copy constructor.
    \param other Another polyhedron
    \return Copied polyhedron
    */
    Polyhedron& operator=(const Polyhedron& other);
    /**
    Polyhedron destructor (deallocates memory of solvers).
    */
    ~Polyhedron();
    /**
    Updates hypeplanes \f$A\in\mathbb{R}^{m\times n}\f$ and \f$b\in\mathbb{R}^{m}\f$. If the number of constraints changes then it deallocates memory for allocated solvers. If the number of contraints are the same. It simply updates the internal variables A_ptr and b_ptr that solvers use.
    \param A Contraints matrix \f$A\in\mathbb{R}^{m\times n}\f$, being \f$m\f$ the number of constraints and \f$n\f$ the polyhedron dimension. Rows of \f$A\f$ are normalized
    \param b Contraints intercepts \f$b\in\mathbb{R}^{m}\f$, being \f$m\f$ the number of constraint
    */
    void update(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
    /**
    Prints the polyhedron constraints on the console.
    */
    virtual void print();
    /**
    Prints the polyhedron constraints on the specified output stream (using Matlab's syntax).
    \param out Output stream to print the polyhedron
    \param AName Name of the variable to assign the constraints matrix \f$A\f$
    \param bName Name of the variable to assign the constraints intercepts \f$b\f$
    \param addRangeLimits If true, it adds to the printed polyhedron the range limits is is always a closed polyhedron
    \return Output stream
    */
    std::ostream& print(std::ostream& out, const std::string &AName, const std::string &bName, const bool& addRangeLimits);
    /**
    Computes constraints from a set of vertices.
    \param[in] v Matrix with \f$v\in\mathbb{R}^{n\times p}\f$, being \f$n\f$ the polyhedron dimension and \f$p\f$ the number of vertices 
    \param[out] A_out Computed contraints matrix \f$A\in\mathbb{R}^{m\times n}\f$, being \f$m\f$ the number of constraints and \f$n\f$ the polyhedron dimension. Rows of \f$A\f$ are normalized
    \param[out] b_out Computed contraints intercepts \f$b\in\mathbb{R}^{m}\f$, being \f$m\f$ the number of constraints
    */
    static void vert2con(const Eigen::MatrixXd& v, Eigen::MatrixXd& A_out, Eigen::VectorXd& b_out);
    /**
    Remove repeated constraints of this polyhedron. If any constraint is removed, the update method is internally executed so that allocated solvers do not fail on next calls.
    \return List with the indices of the removed constraints
    */
    std::vector<int> removeConstraints();
    /**
    Given a constraint matrix and intercepts \f$A\in\mathbb{R}^{m\times n}\f$ and \f$b\in\mathbb{R}^{m}\f$, constraint removal checks if the feasibility of the problem:
    \f{eqnarray*}\displaystyle\min_{x} &\sum x\\ A_{1,i-1}x&\leq b_{1,i-1}\\ A_{i}x&\geq b_{i} \\ A_{i+1,m}x&\leq b_{i+1,m}\f}
    If the problem is feasible then the constraint is redundant and can be removed. It returns the reduced constraint matrices and intercepts \f$A\in\mathbb{R}^{l\times n}\f$ and \f$b\in\mathbb{R}^{l}\f$, with \f$l\leq m\f$.
    \param[in] A Constraints matrix \f$A\in\mathbb{R}^{m\times n}\f$
    \param[in] b Constraints intercepts \f$b\in\mathbb{R}^{m}\f$
    \param[out] Aout Reduced constraints matrix \f$A\in\mathbb{R}^{l\times n}\f$
    \param[out] bout Reduced constraints intercepts \f$b\in\mathbb{R}^{l}\f$
    \return List with the indices of the removed constraints
    */
    static std::vector<int> removeConstraints(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& Aout, Eigen::VectorXd& bout);
    /**
    Compute vertices of a polyhedron from its hyperplane constraints.
    \param[in] v Matrix with \f$v\in\mathbb{R}^{n\times p}\f$, being \f$n\f$ the polyhedron dimension and \f$p\f$ the number of vertices
    \param[out] A Computed contraints matrix \f$A\in\mathbb{R}^{m\times n}\f$, being \f$m\f$ the number of constraints and \f$n\f$ the polyhedron dimension
    \param[out] b Computed contraints intercepts \f$b\in\mathbb{R}^{m}\f$, being \f$m\f$ the number of constraints
    */
    static void con2vert(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::MatrixXd& v);
    //double closestVertex(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    /**
    Computes the closest point \f$x^{*}\f$ of the polyhedron to a given input point \f$x_{in}\f$. The first time this function is called, the solver is allocated and can be used for cached computations for further queries.
    \f{eqnarray*}x^{*}=\text{arg}\!\displaystyle\min_{x} &\|x_{in}-x\|\\ Ax&\leq b\f}
    \param[in] p_in Input point
    \param[out] p_out Computed closest point
    \return Distance of the closest point
    */
    virtual double closestPoint(const Eigen::VectorXd& p_in, Eigen::VectorXd& p_out);
    /**
    Computes the distance between \f$x_{in}\f$ and the closest point \f$x^{*}\f$ of the polyhedron to the given input point. The first time this function is called, the solver is allocated and can be used for cached computations for further queries.
    \f{eqnarray*}x^{*}=\text{arg}\!\displaystyle\min_{x} &\|x_{in}-x\|\\ Ax&\leq b\f}
    \param p_in Input point
    \return Distance of the closest point
    */
    virtual double closestPoint(const Eigen::VectorXd& p_in);
    /**
    Computes the closest point \f$x^{*}\f$ of the polyhedron to a given ellipsoid \f$\mathcal{E}(C,d)=\{x=C\tilde{x}+d\ |\ \|\tilde{x}\|\leq 1\}\f$. The first time this function is called, the solver is allocated and can be used for cached computations for further queries.
    \f{eqnarray*}x^{*}=\text{arg}\!\displaystyle\min_{x,\alpha} &\alpha\\ Ax&\leq b\\ \left[\begin{array}{c}\alpha\\ C^{-1}(x-d)\end{array}\right]\in&\mathcal{Q}^{n+1}\f}
    being \f$\mathcal{Q}^{n}\f$ the domain of quadratic cones \f$\mathcal{Q}^{n}=\{x_1\geq\sqrt{x_2^2+\ldots+x_n^2}\}\f$.
    \param[in] ellipsoid Input ellipsoid \f$\mathcal{E}(C,d)\f$
    \param[out] p_out Computed closest point
    \return Distance of the closest point
    */
    virtual double closestPointExpandingEllipsoid(Ellipsoid& ellipsoid, Eigen::VectorXd& p_out);
    /**
    Check is a point is inside the polyhedron by ensuring that all constraints are satisfied within a given tolerance \f$\delta\f$. Rows of \f$A\f$ are normalized.
    \f[Ax-b\leq \delta \f]
    \param p Point to check if it is inside the polyhedron
    \param tol Tolerance \f$\delta\f$
    \return True if the point is inside the polyhedron, false otherwise.
    */
    virtual bool isInside(const Eigen::VectorXd& p, const double &tol= 1.e-5);
    /**
    Check if part of the polyhedron is inside a given separating hyperplane (defined by the normal vector \f$\tilde{a}_{i}\f$ and intercept \f$\tilde{b}_{i}\f$ and considering a tolerance \f$\delta\f$). This is equivalent to check if the following problem is feasible:
    \f{eqnarray*}\displaystyle\min &\sum_{x} x\\ Ax-b&\leq 0\\ \tilde{a}_{i}x-\tilde{b}_{i}&\leq \delta\f}
    \param ai Normal vector of the separating hyperplane
    \param bi Intercept of the separating hyperplane
    \param tol Tolerance \f$\delta\f$
    \return True if part of the polyhedron is inside the given separating hyperplane
    */
    virtual bool isInsideSeparatingHyperplane(const Eigen::VectorXd& ai, const double& bi, const double &tol= 1.e-5);
    /**
    Check if part of the polyhedron is inside a given set of separating hyperplanes (considering a tolerance \f$\delta\f$). Only the last constraint (last row of \f$A\f$ and \f$b\f$) is indeed considered. The function is cached so that it assumes that will be called every time a new separating hyperplane is added to the set. To reset the solver, you must call this function with an empty constraint set.
    This is equivalent to check if the following problem is feasible:
    \f{eqnarray*}\displaystyle\min &\sum_{x} x\\ Ax-b&\leq 0\\ \tilde{A}x-\tilde{b}&\leq I\cdot\delta\f}
    \param A Constraints matrix of the separating hyperplanes
    \param b Constraints intercepts of the separating hyperplanes
    \param tol Tolerance \f$\delta\f$
    \return True if part of the polyhedron is inside the given separating hyperplane
    */
    virtual bool isInsideSeparatingHyperplanes(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const double& tol = 1.e-5);
    /**
    Returns the Ellipsoid \f$\mathcal{E}(C,c)=\{x=C\tilde{x}+c\ |\ \|\tilde{x}\|\leq 1\}\f$ with the maximum volume inscribed in the polyhedron. The solver is cached and subsequent calls will return the computed ellipsoid if the polyhedron is not changed. This is equivalent to solve the following optimitization problem: 
    \f{eqnarray*}x^{*}=\text{arg}\!\displaystyle\max_{C,c} &\log|C|\\ \left[\begin{array}{c}b-A\cdot c\ A\cdot C\end{array}\right]^{T}\in&\mathcal{Q}^{n+1}\times\mathcal{Q}^{n+1}\ldots\times \mathcal{Q}^{n+1}\\ \left[\begin{array}{c}b_{\mathcal{C}}-A_{\mathcal{C}}\cdot c\ A_{\mathcal{C}}\cdot C\end{array}\right]^{T}\in&\mathcal{Q}^{n+1}\times\mathcal{Q}^{n+1}\ldots\times \mathcal{Q}^{n+1}\f}
    being \f$\mathcal{Q}^{n}\f$ the domain of quadratic cones \f$\mathcal{Q}^{n}=\{x_1\geq\sqrt{x_2^2+\ldots+x_n^2}\}\f$,
    with \f$C\succeq 0\f$ and \f$A_{\mathcal{C}}\in\mathbb{R}^{2n\times n}\f$ and \f$b_{\mathcal{C}}\in\mathbb{R}^{2n}\f$ the separating hyperplanes of the configuration space.
    \return Ellipsoid \f$\mathcal{E}(C,d)=\{x=C\tilde{x}+d\ |\ \|\tilde{x}\|\leq 1\}\f$ inscribed in the polyhedron
    */
    virtual Ellipsoid inscribedEllipsoid();
    /**
    Returns the volume of the Ellipsoid \f$\mathcal{E}(C,c)=\{x=C\tilde{x}+c\ |\ \|\tilde{x}\|\leq 1\}\f$ with the maximum volume inscribed in the polyhedron. The function is a static method, and thus the solver is NOT cached and subsequent calls will allocate and deallocate memory used by the solver. This is equivalent to solve the following optimitization problem:
    \f{eqnarray*}x^{*}=\text{arg}\!\displaystyle\max_{C,c} &\log|C|\\ \left[\begin{array}{c}b-A\cdot c\ A\cdot C\end{array}\right]^{T}\in&\mathcal{Q}^{n+1}\times\mathcal{Q}^{n+1}\ldots\times \mathcal{Q}^{n+1}\\ \left[\begin{array}{c}b_{\mathcal{C}}-A_{\mathcal{C}}\cdot c\ A_{\mathcal{C}}\cdot C\end{array}\right]^{T}\in&\mathcal{Q}^{n+1}\times\mathcal{Q}^{n+1}\ldots\times \mathcal{Q}^{n+1}\f}
    being \f$\mathcal{Q}^{n}\f$ the domain of quadratic cones \f$\mathcal{Q}^{n}=\{x_1\geq\sqrt{x_2^2+\ldots+x_n^2}\}\f$,
    with \f$C\succeq 0\f$ and \f$A_{\mathcal{C}}\in\mathbb{R}^{2n\times n}\f$ and \f$b_{\mathcal{C}}\in\mathbb{R}^{2n}\f$ the separating hyperplanes of the configuration space.
    \param A Constraints matrix of the polyhedron
    \param b Constraints intercepts of the polyhedron
    \return Ellipsoid \f$\mathcal{E}(C,d)=\{x=C\tilde{x}+d\ |\ \|\tilde{x}\|\leq 1\}\f$ inscribed in the polyhedron
    */
    static double inscribedEllipsoidVolume(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
    /**
    Returns the Sphere \f$\mathcal{S}(r,c)=\{x=r\tilde{x}+c\ |\ \|\tilde{x}\|\leq 1\}\f$ with the maximum volume inscribed in the polyhedron. The solver is cached and subsequent calls will return the computed circle if the polyhedron is not changed. This is equivalent to solve the following optimitization problem:
    \f{eqnarray*}x^{*}=\text{arg}\!\displaystyle\max_{r,c} &r\\ \left[\begin{array}{c}b-A\cdot c\ r\cdot A\end{array}\right]^{T}\in&\mathcal{Q}^{n+1}\times\mathcal{Q}^{n+1}\ldots\times \mathcal{Q}^{n+1}\\ \left[\begin{array}{c}b_{\mathcal{C}}-A_{\mathcal{C}}\cdot c\ r\cdot A_{\mathcal{C}}\end{array}\right]^{T}\in&\mathcal{Q}^{n+1}\times\mathcal{Q}^{n+1}\ldots\times \mathcal{Q}^{n+1}\f}
    being \f$\mathcal{Q}^{n}\f$ the domain of quadratic cones \f$\mathcal{Q}^{n}=\{x_1\geq\sqrt{x_2^2+\ldots+x_n^2}\}\f$,
    with \f$r\f$ being the sphere radius and \f$A_{\mathcal{C}}\in\mathbb{R}^{2n\times n}\f$ and \f$b_{\mathcal{C}}\in\mathbb{R}^{2n}\f$ the separating hyperplanes of the configuration space.
    \return Sphere \f$\mathcal{S}(r,c)=\{x=r\tilde{x}+c\ |\ \|\tilde{x}\|\leq 1\}\f$ inscribed in the polyhedron
    */
    virtual Sphere inscribedSphere();
    /**
    Returns the centroid of the inscribed Sphere with the maximum volume in the polyhedron. The computations are cached, meaning that the first time, the function will allocate memory for a solver to compute the sphere, but on subsequent calls, the computed centroid will be returned.
    \return Vector \f$c\in\mathbb{R}^{n}\f$ with the polyhedron centroid
    */
    virtual Eigen::VectorXd getCentroid();
    //virtual bool hasVertices();
    /**
    Checks if two polyhedron intersect. Given \f$\mathcal{P}_{1}(A_{1},b_{1})\f$ and \f$\mathcal{P}_{2}(A_{2},b_{2})\f$ intersect, that is, checks if \f$\left[\begin{array}{c}A_{1}\\A_{2}\end{array}\right]\leq\left[\begin{array}{c}b_{1}\\b_{2}\end{array}\right]\f$ is valid for at least one point. 
    \param polyhedron1 Polyhedron \f$\mathcal{P}_{1}(A_{1},b_{1})\f$
    \param polyhedron2 Polyhedron \f$\mathcal{P}_{2}(A_{2},b_{2})\f$
    \return True if polyhedra intersect
    */
    static bool intersect(const Polyhedron& polyhedron1, const Polyhedron& polyhedron2);
    /**
    Returns the intersection of two polyhedron. Given \f$\mathcal{P}_{1}(A_{1},b_{1})\f$ and \f$\mathcal{P}_{2}(A_{2},b_{2})\f$ return \f$\left[\begin{array}{c}A_{1}\\A_{2}\end{array}\right]\leq\left[\begin{array}{c}b_{1}\\b_{2}\end{array}\right]\f$.
    \param polyhedron1 Polyhedron \f$\mathcal{P}_{1}(A_{1},b_{1})\f$
    \param polyhedron2 Polyhedron \f$\mathcal{P}_{2}(A_{2},b_{2})\f$
    \return Polyhedron intersection
    */
    static Polyhedron intersection(const Polyhedron& polyhedron1, const Polyhedron& polyhedron2);
    /**
    Checks if the current polyhedron intersects with another polyhedron. Given another polyhedron \f$\mathcal{P}'(A',b')\f$, checks if \f$\left[\begin{array}{c}A\\A'\end{array}\right]\leq\left[\begin{array}{c}b\\b'\end{array}\right]\f$ is valid for at least one point.
    \param A Constraint matrix of the (another) polyhedron
    \param b Contraint intercept of the (another) polyhedron
    \return True if polyhedron intersects with the current polyhedron
    */
    bool intersect(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
    /**
    Returns the bounding box of the polyhedron. The function is cached, meaning that the first time is called computes the bounding box, while subsequent calls will return the computed bounding box if the polyhedron is not changed. The bounding box is defined as a set of vertices.
    \return Bounding box vertices with dimensions \f$n\times 2n\f$, being \f$n\f$ the dimension of the polyhedron 
    */
    virtual Eigen::MatrixXd getBoundingBox();
    /**
    Returns the constraint that is closest to the point \f$p\in\mathcal{P}(A,b)\f$.
    \param[in] p Point belonging to the polyhedron
    \param[out] ai Constraint normal vector of the constraint closest to the point
    \param[out] bi Constraint intercept of the constraint closest to the point 
    */
    void closestConstraint(const Eigen::VectorXd& p, Eigen::VectorXd& ai, double& bi);
    /**
    Returns the indexes of constraints \f$\{a_{i},b_{i}\}\f$ that satisfy \f$|a_{i}p-b_{i}|\delta\f$ for a given point belongin to the polyhedron boundary with tolerance \f$\delta\f$.
    \param p Point \f$p\f$ belonging to the polyhedron boundary
    \param tol Tolerance \f$\delta\f$
    \return Vector with the indexes of the contraints satisfying the condition 
    */
    std::vector<int> eqConstraints(const Eigen::VectorXd& p, const double& tol = 1.e-5);
    void deallocateInscribedEllipsoidSolver();
protected:
    //Polyhedron constraints
    std::shared_ptr<ndarray<double, 2>> A_ptr;
    std::shared_ptr<ndarray<double, 1>> b_ptr;
    Model::t MClosestPoint;
    Variable::t xClosestPoint;
    Variable::t l;
    Ellipsoid ellipsoid;
    bool computedInscribedEllipsoid=false;
    Sphere sphere;
    bool computedInscribedSphere=false;
    bool solverClosestPointAllocated =false;
    Model::t MClosestPointEllipsoid;
    Variable::t xClosestPointEllipsoid;
    Variable::t wClosestPointEllipsoid;
    Variable::t alphaClosestPointEllipsoid;
    bool solverClosestPointEllipsoidAllocated = false;
    Model::t MIsInsideSeparatingHyperplane;
    Variable::t xIsInsideSeparatingHyperplane;
    bool solverIsInsideSeparatingHyperplaneAllocated = false;
    int numIsInsideSeparatingHyperplanes = 0;
    void allocateClosestPointEllipsoidQPSolver();
    Model::t MInscribedEllipsoid;
    Variable::t xInscribedEllipsoid;
    Variable::t CInscribedEllipsoid;
    Variable::t tInscribedEllipsoid;
    Variable::t YInscribedEllipsoid;
    Variable::t ZInscribedEllipsoid;
    Variable::t DZInscribedEllipsoid;
    bool solverInscribedEllipsoidAllocated = false;
    void computeInscribedEllipsoid();
    void computeInscribedSphere();
    virtual void allocateClosestPointSolver();
    virtual void allocateClosestPointEllipsoidSolver();
    virtual void allocateIsInsideSeparatingHyperplaneSolver();
    void allocateInscribedEllipsoidSolver();
};
#endif