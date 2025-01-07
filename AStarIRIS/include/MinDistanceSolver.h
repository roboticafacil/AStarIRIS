#ifndef MIN_DISTANCE_SOLVER
#define MIN_DISTANCE_SOLVER
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
}MIPSol_t;

class MinDistanceSolver : public Solver
{
public:
	MIPSol_t MIPSolution;
	MinDistanceSolver(const Graph& g, const int& startKey, const int& targetKey);
	void setTask();
private:
	Variable::t l;
	Variable::t z;
	Variable::t p;
	Variable::t y;
	//Variable::t v;
public:
	void solve();
private:
	void computeFeasibleSolution(const int &maxIters=1);
};
#endif