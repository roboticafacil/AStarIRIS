#ifndef CONVEX_RELAXATION
#define CONVEX_RELAXATION
#include "Graph.h"
#include "fusion.h"
#include "EigenNdArray.h"

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

typedef struct
{
	double cost;
	Eigen::MatrixXd x;
	std::vector<int> nodeKeys;
	std::vector<int> edgeKeys;
}FeasibleSol_t;

class ConvexRelaxation
{
public:
	ConvexRelaxationSol_t relaxedSolution;
	FeasibleSol_t feasibleSolution;
protected:
	Graph g;
	Model::t M;
	Variable::t l;
	Variable::t z;
	Variable::t p;
	Variable::t y;
	Variable::t v;
	int n;
	int nE;
	int nN;
	int startKey;
	int targetKey;
	std::shared_ptr<ndarray<double, 2>> qstart;
	std::shared_ptr<ndarray<double, 2>> qtarget;
	
public:
	ConvexRelaxation(const Graph& g, const int& startKey, const int& targetKey);
	virtual void setTask()=0;
	void solve();
private:
	void computeFeasibleSolution();
};
#endif