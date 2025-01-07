#ifndef SOLVER
#define SOLVER
#include "Graph.h"
#include "fusion.h"
#include "EigenNdArray.h"

using namespace mosek::fusion;
using namespace monty;

typedef struct
{
	std::vector<int> nodeKeys;
	std::vector<int> edgeKeys;
}Path_t;

typedef struct
{
	SolutionStatus status;
	double cost;
	Eigen::MatrixXd x;
}FeasibleSol_t;

class Solver
{
public:
	FeasibleSol_t feasibleSolution;
	Path_t optimalPath;
protected:
	Graph g;
	Model::t M;
	int n;
	int nE;
	int nN;
	int startKey;
	int targetKey;
	std::shared_ptr<ndarray<double, 2>> qstart;
	std::shared_ptr<ndarray<double, 2>> qtarget;

public:
	Solver(const Graph& g, const int& startKey, const int& targetKey);
	~Solver();
	virtual void setTask() = 0;
	virtual void solve() = 0;
	void setGraph(Graph& g);
	bool optimalPathContainsTerminalNodes();
protected:
	virtual void computeFeasibleSolution(const int &maxIters=1)=0;
	bool solverAllocated;
};

#endif