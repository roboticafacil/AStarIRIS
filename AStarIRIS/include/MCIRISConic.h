#ifndef MC_IRIS_CONIC
#define MC_IRIS_CONIC
#include "Graph.h"
#include<Eigen/Dense>
#include "CObsConic.h"
#include "GCS.h"
#include "ExpandableIRISConic.h"
#include "PointNode.h"
#include "PerspectiveSPP_GCS.h"

typedef struct
{
	ExpandableIRISParams_t ExpandableIRISParams;
	//Add here more parameters if needed
	int maxIters;
}MCIRISParams_t;

class MCIRISConic : public ExpandableIRISConic
{
public:
	FeasibleSol_t feasibleSolution; //These two properties probably should be inherited from a base class (see AStarIRISConic also has the same properties)
	Path_t optimalPath;
public:
	MCIRISConic(CObsConic& cObs, const MCIRISParams_t& params);
	static MCIRISParams_t getDefaultMCIRISParams();
	void do_MCRelaxedSolver();
	MCIRISParams_t params;
};
#endif