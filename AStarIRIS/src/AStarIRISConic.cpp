#include "AStarIRISConic.h"
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


AStarIRISConic::AStarIRISConic(CObsConic& cObs, const AStarIRISParams_t& params) : ExpandableIRISConic(cObs,params.ExpandableIRISParams), params(params)
{

}


void AStarIRISConic::doPhase1_RelaxedSolver()
{
	this->addConvexSets(qStartNode->point.p);
	bool phase1_finished = this->gcs.contains(qTargetNode->point.p);
	while (!phase1_finished)
	{
		ConvexRelaxationMinDistanceSolver solver(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solver.setTask();
		solver.solve();
		//solver.getGreedyPath();
		solver.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);
		//Select the terminal node connected to the target that generated the feasible solution (the size of nodeKeys is always at least 2, because it includes the start and target nodes)
		int terminalNodeKey=solver.optimalPath.nodeKeys[solver.optimalPath.nodeKeys.size()-2];
		Node* node = this->navGraph.getNode(terminalNodeKey);
		if (dynamic_cast<PolyhedronTerminalNode*>(node) == NULL)
		{
			std::cout << "We have found a feasible solution" << std::endl;
			//solver.computeFeasibleSolution();
			//this->feasibleSolution = solver.feasibleSolution;
			break;
		}
		std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;
		this->expandTerminalNode(terminalNodeKey);
		/*std::cout << "Graph of convex sets" << std::endl;
		std::vector<int> nodeKeys = this->gcs.getNodeKeys();
		for (int i = 0; i < this->gcs.numNodes; i++)
		{
			std::cout << "Convex set " << nodeKeys[i] << std::endl;
			Node* node = this->gcs.getNode(nodeKeys[i]);
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		this->gcs.print();
		std::cout << "NavGraph of convex sets" << std::endl;
		std::vector<int> edgeNodeKeys = this->navGraph.getNodeKeys();
		for (int i = 0; i < this->navGraph.numNodes; i++)
		{
			std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
			Node* node = this->navGraph.getNode(edgeNodeKeys[i]);
			if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
				PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
				polyNode->polyhedron.print();
			}
			else if (dynamic_cast<PointNode*>(node) != NULL) {
				PointNode* pointNode = (PointNode*)node->getNodeData();
				pointNode->point.print();
			}
		}
		this->navGraph.print();
		std::cout << "Relaxed Solution " << std::endl;
		std::cout << "x variables " << std::endl;
		std::cout << solver.relaxedSolution.x << std::endl;
		std::cout << "y variables " << std::endl;
		std::cout << solver.relaxedSolution.y << std::endl;
		std::cout << "z variables " << std::endl;
		std::cout << solver.relaxedSolution.z << std::endl;
		std::cout << "p variables " << std::endl;
		std::cout << solver.relaxedSolution.p << std::endl;
		std::cout << "l variables " << std::endl;
		std::cout << solver.relaxedSolution.l << std::endl;
		std::cout << "cost " << std::endl;
		std::cout << solver.relaxedSolution.cost << std::endl;*/
		this->optimalPath = solver.optimalPath;
		/*std::cout << "Feasible Solution " << std::endl;
		std::cout << "x variables " << std::endl;
		std::cout << solver.feasibleSolution.x << std::endl;
		std::cout << "nodes keys ";
		for (std::vector<int>::iterator it = solver.feasibleSolution.nodeKeys.begin(); it != solver.feasibleSolution.nodeKeys.end(); it++)
			std::cout << *it << " ";
		std::cout << std::endl;
		std::cout << "edge keys ";
		for (std::vector<int>::iterator it = solver.feasibleSolution.edgeKeys.begin(); it != solver.feasibleSolution.edgeKeys.end(); it++)
			std::cout << *it << " ";
		std::cout << std::endl;
		std::cout << "cost " << std::endl;
		std::cout << solver.feasibleSolution.cost << std::endl;*/
		/*phase1_finished = true;
		std::vector<int> nodeKeys = this->navGraph.getNodeKeys();
		for (std::vector<int>::iterator it = nodeKeys.begin(); it != nodeKeys.end(); it++)
		{
			Node* node = (Node*)this->navGraph.getNode(*it);
			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				phase1_finished = false;
				break;
			}
		}*/
		phase1_finished = this->gcs.contains(qTargetNode->point.p);
		//if (phase1_finished)
		//{
		//	solver.computeFeasibleSolution();
		//	this->feasibleSolution = solver.feasibleSolution;
		//}
	}
}

void AStarIRISConic::doPhase1_MIPSolver()
{
	this->addConvexSets(qStartNode->point.p);
	bool phase1_finished = this->gcs.contains(qTargetNode->point.p);
	while (!phase1_finished)
	{
		MinDistanceSolver solver(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solver.setTask();
		solver.solve();
		//Select the terminal node connected to the target that generated the feasible solution (the size of nodeKeys is always at least 2, because it includes the start and target nodes)
		int terminalNodeKey = solver.optimalPath.nodeKeys[solver.optimalPath.nodeKeys.size() - 2];
		std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;
		//if (terminalNodeKey == 38)
		//{
		//	std::cout << "This is the one" << std::endl; //For some strange reason the seed generator is not working properly, even if I increase the number of maximum trials, but apparently the facet is expandable
		//	//params.maxTrialsTerminal = 999;
		//}
		this->expandTerminalNode(terminalNodeKey);
		

		/*std::vector<int> nodeKeysFound = this->gcs.findConvexSets(qTargetNode->point.p);
		bool qTargetInGCS = nodeKeysFound.size()>0;
		if (qTargetInGCS)
		{
			std::cout << "Graph of convex sets" << std::endl;
			std::vector<int> nodeKeys = this->gcs.getNodeKeys();
			for (int i = 0; i < this->gcs.numNodes; i++)
			{
				std::cout << "Convex set " << nodeKeys[i] << std::endl;
				Node* node = this->gcs.getNode(nodeKeys[i]);
				PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
				polyNode->polyhedron.print();
			}
			this->gcs.print();
			std::cout << "NavGraph of convex sets" << std::endl;
			std::vector<int> edgeNodeKeys = this->navGraph.getNodeKeys();
			for (int i = 0; i < this->navGraph.numNodes; i++)
			{
				Node* node = this->navGraph.getNode(edgeNodeKeys[i]);
				if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL) {
					std::cout << "Terminal Node Convex set " << edgeNodeKeys[i] << std::endl;
					PolyhedronTerminalNode* polyNode = (PolyhedronTerminalNode*)node->getNodeData();
					polyNode->polyhedron.print();
				}
				else if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
					std::cout << "Polyhedron Node Convex set " << edgeNodeKeys[i] << std::endl;
					PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
					polyNode->polyhedron.print();
				}
				else if (dynamic_cast<PointNode*>(node) != NULL) {
					std::cout << "Point Node Convex set " << edgeNodeKeys[i] << std::endl;
					PointNode* pointNode = (PointNode*)node->getNodeData();
					pointNode->point.print();
				}
			}
			this->navGraph.print();
			std::cout << "MIP Solution " << std::endl;
			std::cout << "x variables " << std::endl;
			std::cout << solver.MIPSolution.x << std::endl;
			std::cout << "y variables " << std::endl;
			std::cout << solver.MIPSolution.y << std::endl;
			std::cout << "z variables " << std::endl;
			std::cout << solver.MIPSolution.z << std::endl;
			std::cout << "p variables " << std::endl;
			std::cout << solver.MIPSolution.p << std::endl;
			std::cout << "l variables " << std::endl;
			std::cout << solver.MIPSolution.l << std::endl;
			std::cout << "cost " << std::endl;
			std::cout << solver.MIPSolution.cost << std::endl;
			this->feasibleSolution = solver.feasibleSolution;
			std::cout << "Feasible Solution " << std::endl;
			std::cout << "x variables " << std::endl;
			std::cout << solver.feasibleSolution.x << std::endl;
			std::cout << "nodes keys ";
			for (std::vector<int>::iterator it = solver.feasibleSolution.nodeKeys.begin(); it != solver.feasibleSolution.nodeKeys.end(); it++)
				std::cout << *it << " ";
			std::cout << std::endl;
			std::cout << "edge keys ";
			for (std::vector<int>::iterator it = solver.feasibleSolution.edgeKeys.begin(); it != solver.feasibleSolution.edgeKeys.end(); it++)
				std::cout << *it << " ";
			std::cout << std::endl;
			std::cout << "cost " << std::endl;
			std::cout << solver.feasibleSolution.cost << std::endl;
		}
		Eigen::VectorXd previousConf=solver.feasibleSolution.x.row(solver.feasibleSolution.x.rows()-2);
		bool previousConfValid = false;
		for (std::vector<int>::iterator itNodeKeysFound = nodeKeysFound.begin(); itNodeKeysFound != nodeKeysFound.end(); itNodeKeysFound++)
		{
			Node* node=(Node*)this->gcs.getNode(*itNodeKeysFound);
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			if (polyNode->polyhedron.isInside(previousConf))
			{
				previousConfValid = true;
				break;
			}
		}
		phase1_finished = qTargetInGCS && previousConfValid;*/
		/*phase1_finished = true;
		std::vector<int> nodeKeys = this->navGraph.getNodeKeys();
		for (std::vector<int>::iterator it = nodeKeys.begin(); it != nodeKeys.end(); it++)
		{
			Node* node = (Node*)this->navGraph.getNode(*it);
			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				phase1_finished = false;
				break;
			}
		}*/
		phase1_finished = this->gcs.contains(qTargetNode->point.p);
	}
}


AStarIRISParams_t AStarIRISConic::getDefaultAStarIRISParams()
{
	AStarIRISParams_t params = { ExpandableIRISConic::getDefaultExpandableIRISParams()};
	return params;
}