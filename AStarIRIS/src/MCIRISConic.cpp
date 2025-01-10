#include "MCIRISConic.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "EllipsoidNode.h"
#include "IRISConic.h"
#include "PolyhedronObstacleCircularRobotNode.h"
#include "CObsConic.h"
#include "GCS.h"
#include "ConvexRelaxationMinDistanceSolver.h"
#include "MinDistanceSolver.h"


MCIRISConic::MCIRISConic(CObsConic& cObs, const MCIRISParams_t& params) : ExpandableIRISConic(cObs, params.ExpandableIRISParams), params(params)
{

}

void MCIRISConic::do_MCRelaxedSolver()
{
	this->addConvexSets(qStartNode->point.p);
	bool targetFound = this->gcs.contains(qTargetNode->point.p);
	int iters = 0;
	double bestCost=(double)std::numeric_limits<double>::infinity();
	//Phase 1
	while (!this->gcs.contains(qTargetNode->point.p))
	{
		ConvexRelaxationMinDistanceSolver solver(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solver.setTask();
		solver.solve();
		solver.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);
		int previousLastNodeKey = solver.optimalPath.nodeKeys[solver.optimalPath.nodeKeys.size() - 2];
		Node* node = this->navGraph.getNode(previousLastNodeKey);
		if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
		{
			std::cout << "Terminal node to expand " << previousLastNodeKey << std::endl;
			this->expandTerminalNode(previousLastNodeKey);
			targetFound = this->gcs.contains(qTargetNode->point.p);
		}
		else
		{
			std::cout << "Warning:: This shouldn't happen, because we haven't found the GCS node that contains the target, so all connections to target must be with terminal nodes" << std::endl;
		}
	}
	//Phase 2
	NavGraph simplifiedGraphPhase1 = this->getGraphWithoutTerminalConnections();
	ConvexRelaxationMinDistanceSolver feasibleSolverPhase1(simplifiedGraphPhase1, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
	feasibleSolverPhase1.setTask();
	feasibleSolverPhase1.solve();
	feasibleSolverPhase1.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);
	bestCost = feasibleSolverPhase1.feasibleSolution.cost;
	this->optimalPath = feasibleSolverPhase1.optimalPath;
	this->feasibleSolution = feasibleSolverPhase1.feasibleSolution;
	while (iters < params.maxIters)
	{
		ConvexRelaxationMinDistanceSolver solver(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solver.setTask();
		solver.solve();
		solver.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);
		if (solver.feasibleSolution.cost < bestCost)
		{
			int previousLastNodeKey = solver.optimalPath.nodeKeys[solver.optimalPath.nodeKeys.size() - 2];
			Node* node = this->navGraph.getNode(previousLastNodeKey);
			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				std::cout << "Terminal node to expand " << previousLastNodeKey << std::endl;
				this->expandTerminalNode(previousLastNodeKey);
			}
			NavGraph simplifiedGraph = this->getGraphWithoutTerminalConnections();
			ConvexRelaxationMinDistanceSolver feasibleSolver(simplifiedGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
			feasibleSolver.setTask();
			feasibleSolver.solve();
			feasibleSolver.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);
			if (feasibleSolver.feasibleSolution.cost < bestCost)
			{
				bestCost = feasibleSolver.feasibleSolution.cost;
				this->optimalPath = feasibleSolver.optimalPath;
				this->feasibleSolution = feasibleSolver.feasibleSolution;
				iters = 0;
			}
		}
		else
		{
			break;
		}
		iters++;
	}
}


MCIRISParams_t MCIRISConic::getDefaultMCIRISParams()
{
	MCIRISParams_t params = { ExpandableIRISConic::getDefaultExpandableIRISParams(),15};
	return params;
}