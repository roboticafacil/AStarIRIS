#ifndef CONVEX_RELAXATION_MIN_DISTANCE
#define CONVEX_RELAXATION_MIN_DISTNACE
#include "Graph.h"
#include "fusion.h"
#include "EigenNdArray.h"
#include "ConvexRelaxation.h"

using namespace mosek::fusion;
using namespace monty;

class ConvexRelaxationMinDistance: public ConvexRelaxation
{
public:
	ConvexRelaxationMinDistance(const Graph& g, const int& startKey, const int& targetKey);
	void setTask();
};
#endif