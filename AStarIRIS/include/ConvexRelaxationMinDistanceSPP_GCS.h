#ifndef CONVEX_RELAXATION_MIN_DISTANCE_SPP_GCS_H
#define CONVEX_RELAXATION_MIN_DISTANCE_SPP_GCS_H
#include "Graph.h"
#include "fusion.h"
#include "EigenUtils.h"
#include "PerspectiveSPP_GCS.h"

using namespace mosek::fusion;
using namespace monty;

class ConvexRelaxationMinDistanceSPP_GCS: public PerspectiveSPP_GCS
{
public:
	ConvexRelaxationMinDistanceSPP_GCS(Graph* g, const int& startKey, const int& targetKey);
	void optimize();
	void computeFeasibleSolution(const int& maxIters=1, const bool& simplifyGraph=false);
protected:
	void allocateVariables(Model::t& M);
	void setConstraints(Model::t& M);
	void setObjective(Model::t& M);
};
#endif