#ifndef A_STAR_IRIS_CONIC
#define A_STAR_IRIS_CONIC
#include "Graph.h"
#include<Eigen/Dense>
#include "CObsConic.h"
#include "GCS.h"
#include "ExpandableIRISConic.h"
#include "PointNode.h"
#include "PerspectiveSPP_GCS.h"
#include "ConvexRelaxationMinDistanceSPP_GCS.h"
#include "ConvexRelaxationBezierCurveSPP_GCS.h"

typedef struct
{
	ExpandableIRISParams_t ExpandableIRISParams;
	//Add here more parameters if needed
}AStarIRISParams_t;

typedef struct
{
	bool gcsGraph = false;
	bool navGraph = false;
	bool gcsConvexSets = false;
	bool navGraphConvexSets = false;
	bool video_frames = false;
	bool phase1 = false;
	bool phase2 = false;
	int multiagent = 1;
}AStartIRISDebugLevel_t;

class AStarIRISConic: public ExpandableIRISConic
{
public:
	FeasibleSol_t feasibleSolution;
	Path_t optimalPath;
public:
	AStarIRISConic(CObsConic& cObs, const AStarIRISParams_t& params);
	static AStarIRISParams_t getDefaultAStarIRISParams();
	void do_RelaxedSolver(PerspectiveSPP_GCS& solver);
	void do_RelaxedSolver(PerspectiveSPP_GCS& solver, std::ostream& out, const AStartIRISDebugLevel_t debug_level);
	void do_MIPSolver();
	void do_MIPSolver(std::ostream& out);
	AStarIRISParams_t params;
};
#endif