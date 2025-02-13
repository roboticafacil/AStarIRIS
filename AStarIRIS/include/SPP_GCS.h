#ifndef SPP_GCS_H
#define SPP_GCS_H
#include "Graph.h"
#include "fusion.h"
#include "EigenUtils.h"

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

class SPP_GCS
{
public:
	FeasibleSol_t feasibleSolution;
	Path_t optimalPath;
protected:
	Graph* g;
	//Model::t M;
	int n;
	//int nE;
	//int nN;
	int N;
	int startKey;
	int targetKey;
	std::shared_ptr<ndarray<double, 2>> qstart;
	std::shared_ptr<ndarray<double, 2>> qtarget;
	//Eigen::MatrixXi idx;

public:
	SPP_GCS(Graph* g, const int &N, const int& startKey, const int& targetKey);
	virtual void optimize() = 0;
	void setGraph(Graph* g);
	bool optimalPathContainsTerminalNodes();
	virtual void computeFeasibleSolution(const int& maxIters = 1, const bool& simplifyGraph = false) = 0;
};

#endif