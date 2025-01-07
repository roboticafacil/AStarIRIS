#ifndef CONVEX_RELAXATION_MIN_DISTANCE_SOLVER
#define CONVEX_RELAXATION_MIN_DISTANCE_SOLVER
#include "Graph.h"
#include "fusion.h"
#include "EigenNdArray.h"
#include "ConvexRelaxationSolver.h"

using namespace mosek::fusion;
using namespace monty;

class ConvexRelaxationMinDistanceSolver: public ConvexRelaxationSolver
{
public:
	ConvexRelaxationMinDistanceSolver(const Graph& g, const int& startKey, const int& targetKey);
	void setTask();
	void computeFeasibleSolution(const int& maxIters=1);
};
#endif