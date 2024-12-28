#ifndef CONVEX_RELAXATION_MIN_DISTANCE
#define CONVEX_RELAXATION_MIN_DISTNACE
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
	Eigen::MatrixXd zp;
	Eigen::MatrixXd x;
}ConvexRelaxationSol_t;

class ConvexRelaxationMinDistance
{
private:
	Graph g;
	std::vector<int> nodeKeys;
	Model::t M;
	Variable::t l;
	Variable::t z;
	Variable::t zp;
	Variable::t y;
	Variable::t v;
	int n;
	int nE;
	int nN;
public:
	ConvexRelaxationMinDistance(const Graph& g, const std::vector<int>& nodeKeys, const int& startKey, const int& targetKey);
	void setMinDistanceTask();
	ConvexRelaxationSol_t solve();
private:
	int startKey;
	int targetKey;
	std::shared_ptr<ndarray<double, 2>> qstart;
	std::shared_ptr<ndarray<double, 2>> qtarget;
};
#endif