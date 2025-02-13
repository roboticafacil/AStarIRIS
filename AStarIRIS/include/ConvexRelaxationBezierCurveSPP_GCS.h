#ifndef CONVEX_RELAXATION_BEZIER_CURVE_SPP_GCS_H
#define CONVEX_RELAXATION_BEIZER_CURVE_SPP_GCS_H
#include "Graph.h"
#include "fusion.h"
#include "EigenUtils.h"
#include "PerspectiveSPP_GCS.h"

using namespace mosek::fusion;
using namespace monty;

typedef struct
{
	int N;
	int costDerivativeOrder;
	int continuityDerivativeOrder;
	bool squaredRootL2NormCost;
	bool useNavGraph;
}ConvexRelaxationBezierCurveParams_t;


class ConvexRelaxationBezierCurveSPP_GCS : public PerspectiveSPP_GCS
{
public:
	ConvexRelaxationBezierCurveSPP_GCS(Graph* navGraph, Graph* gcs, std::map<std::pair<int, int>, int>* navGraph2gcs, const ConvexRelaxationBezierCurveParams_t& params, const int& startKey, const int& targetKey, const Eigen::MatrixXd& CpStart, const Eigen::MatrixXd& CpTarget);
	void optimize();
	void computeFeasibleSolution(const int& maxIters = 1, const bool& simplifyGraph = false);
	static ConvexRelaxationBezierCurveParams_t getDefaultBezierCurveSolverParams();
protected:
	std::shared_ptr<ndarray<double, 1>>  cpStart;
	std::shared_ptr<ndarray<double, 1>>  cpTarget;
	std::shared_ptr<ndarray<double, 1>>  cppStart;
	std::shared_ptr<ndarray<double, 1>>  cppTarget;
	std::shared_ptr<ndarray<double, 1>>  cpppStart;
	std::shared_ptr<ndarray<double, 1>>  cpppTarget;
	Graph* gcs;
	std::map<std::pair<int, int>, int>* navGraph2gcs;
	ConvexRelaxationBezierCurveParams_t params;
	std::shared_ptr<ndarray<double, 2>> Qptr;
	Eigen::MatrixXd computeBezierL2NormMatrix();
	void allocateVariables(Model::t& M);
	void setConstraints(Model::t& M);
	void setObjective(Model::t& M);
};
#endif