#ifndef MIP_MIN_DISTANCE_SPP_GCS_H
#define MIP_MIN_DISTANCE_SPP_GCS_H
#include "Graph.h"
#include "fusion.h"
#include "EigenUtils.h"
#include "ConvexRelaxationMinDistanceSPP_GCS.h"

using namespace mosek::fusion;
using namespace monty;

typedef struct
{
	SolutionStatus status;
	double cost;
	int N;
	Eigen::VectorXd y;
	Eigen::MatrixXd z;
	Eigen::MatrixXd x;
	Eigen::MatrixXd l;
}MIPSol_t;

class MIPMinDistanceSPP_GCS : public ConvexRelaxationMinDistanceSPP_GCS
{
public:
	//MIPSol_t MIPSolution;
	MIPMinDistanceSPP_GCS(Graph* g, const int& startKey, const int& targetKey);
private:
	void allocateVariables(Model::t& M);
};
#endif