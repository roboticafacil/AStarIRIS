#ifndef CONVEX_RELAXATION_SOLVER
#define CONVEX_RELAXATION_SOLVER
#include "Graph.h"
#include "fusion.h"
#include "EigenNdArray.h"
#include "Solver.h"

using namespace mosek::fusion;
using namespace monty;

typedef struct
{
	SolutionStatus status;
	double cost;
	Eigen::VectorXd y;
	Eigen::MatrixXd z;
	Eigen::MatrixXd p;
	Eigen::MatrixXd x;
	Eigen::MatrixXd l;
}ConvexRelaxationSol_t;


class ConvexRelaxationSolver: public Solver
{
public:
	ConvexRelaxationSol_t relaxedSolution;
protected:
	Variable::t l;
	Variable::t z;
	Variable::t p;
	Variable::t y;
	//Variable::t v;
	
public:
	ConvexRelaxationSolver(const Graph& g, const int& startKey, const int& targetKey);
	virtual void setTask()=0;
	void solve();
	Path_t getGreedyPath();
	Path_t getMCPath();
protected:
	virtual void computeFeasibleSolution(const int &maxIters=1) = 0;
};
#endif