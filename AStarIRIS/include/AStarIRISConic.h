#ifndef A_STAR_IRIS_CONIC
#define A_STAR_IRIS_CONIC
#include "Graph.h"
#include<Eigen/Dense>
#include "CObsConic.h"
#include "GCS.h"
#include "ExpandableIRISConic.h"
#include "PointNode.h"
#include "ConvexRelaxationSolver.h"

typedef struct
{
	ExpandableIRISParams_t ExpandableIRISParams;
	//Add here more parameters if needed
}AStarIRISParams_t;

class AStarIRISConic: public ExpandableIRISConic
{
public:
	FeasibleSol_t feasibleSolution;
	Path_t optimalPath;
public:
	AStarIRISConic(CObsConic& cObs, const AStarIRISParams_t& params);
	static AStarIRISParams_t getDefaultAStarIRISParams();
	void doPhase1_RelaxedSolver();
	void doPhase1_MIPSolver();
	AStarIRISParams_t params;
};
#endif