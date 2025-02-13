#ifndef PERSPECTIVE_SPP_GCS
#define PERSPECTIVE_SPP_GCS
#include "Graph.h"
#include "fusion.h"
#include "EigenUtils.h"
#include "SPP_GCS.h"

using namespace mosek::fusion;
using namespace monty;

typedef struct
{
	SolutionStatus status;
	double cost;
	int N;
	Eigen::VectorXd y;
	Eigen::MatrixXd z;
	Eigen::MatrixXd l;
}PerspectiveSol_t;


class PerspectiveSPP_GCS: public SPP_GCS
{
public:
	PerspectiveSol_t perspectiveSolution;
protected:
	Variable::t l;
	Variable::t z;
	Variable::t y;
	
public:
	PerspectiveSPP_GCS(Graph* g, const int &N, const int& startKey, const int& targetKey);
	virtual void optimize()=0;
	Path_t getGreedyPath();
	bool getMCPath(Path_t& path);
	void simplifyGraph();
	virtual void computeFeasibleSolution(const int& maxIters = 1, const bool& simplifyGraph = false) = 0;
	virtual void allocateVariables(Model::t& M)=0;
	virtual void setConstraints(Model::t& M)=0;
	virtual void setObjective(Model::t& M) = 0;
protected:	
	Graph simplifiedGraph;
	bool graphSimplified;
	bool perspectiveSolutionSolved;
};
#endif