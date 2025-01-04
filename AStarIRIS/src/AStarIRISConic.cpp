#include "AStarIRISConic.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "EllipsoidNode.h"
#include "IRISConic.h"
#include "PolyhedronObstacleCircularRobotNode.h"
#include "CObsConic.h"
#include "GCS.h"
#include <ConvexRelaxationMinDistance.h>

/*AStarIRISConic::AStarIRISConic(const Graph& g, const std::vector<int>& nodeKeys, const int& startKey, const int& targetKey) : graph(g), nodeKeys(nodeKeys), h(g.numNodes), nN(g.numNodes), startKey(startKey), targetKey(targetKey)
{
	PointNode* node = (PointNode*)graph.getNode(startKey);
	this->n = node->point.p.size();
	this->pClosest=Eigen::MatrixXd(n,nN);
	this->computeHeuristics();
}*/

AStarIRISConic::AStarIRISConic(CObsConic& cObs, const AStarIRISParams_t& params) : IRISConic(cObs,params.IRISParams), params(params)
{

}

void AStarIRISConic::addConvexSets(const Eigen::VectorXd& q)
{
	Ellipsoid ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	Polyhedron* convexSet = new Polyhedron(q.rows());
	std::vector<IRISNeighbour_t> neighbours;
	//std::vector<int> newNodeKeys;
	//convexSet.allocateClosestPointEllipsoidSolver();
	//convexSet.allocateInscribedEllipsoidSolver();
	while (true)
	{
		IRISConic::computeConvexSet(ellipsoid, *convexSet, neighbours);
		//Add the convexSet to GCS
		Eigen::VectorXd bShrinked = (1. - this->params.IRISParams.shrinkFactor) * convexSet->b + this->params.IRISParams.shrinkFactor * (convexSet->A * ellipsoid.getCentroid());
		PolyhedronNode* node = new PolyhedronNode(convexSet->A, convexSet->b, bShrinked);
		bool isStartInside = convexSet->isInside(qStartNode->point.p);
		//bool isTargetInside = convexSet->isInside(qTargetNode->point.p);
		int nodeKey = this->gcs.addNode(node);
		//newNodeKeys.push_back(nodeKey);
		//Add terminal nodes to EGCS
		for (int i = 0; i < convexSet->A.rows(); i++)
		{
			Eigen::MatrixXd facetA(convexSet->A.rows() + 1, convexSet->A.cols());
			Eigen::VectorXd facetB(convexSet->A.rows() + 1);
			facetA << convexSet->A, -convexSet->A.row(i);
			facetB << convexSet->b, -convexSet->b(i) + 1e-3;
			//int terminalEGCSKey = this->egcs.addNode(new PolyhedronTerminalNode(facetA, facetB,i));
			int terminalEGCSKey = this->navGraph.addNode(new PolyhedronTerminalNode(facetA, facetB, i));
			terminalEGCSNodeKeys.push_back(std::make_pair(nodeKey, terminalEGCSKey));
			if (isStartInside)
				this->navGraph.addEdge(this->qStartNodeNavGraphKey, terminalEGCSKey);
			this->navGraph.addEdge(terminalEGCSKey, this->qTargetNodeNavGraphKey);
		}
		//Add neighbours to GCS, create EGCS node, set the EGCS to GCS map and set neighbours of EGCS node
		std::vector<int> neighbourNodeEGCSKeys;
		for (std::vector<IRISNeighbour_t>::iterator it = neighbours.begin(); it != neighbours.end(); it++)
		{
			//Add neighbours to GCS
			this->gcs.addEdge(nodeKey, it->nodeKey);
			//this->gcs->addEdge(it->nodeKey,nodeKey);  //Maybe this is not necessary
			//Create EGCS node
			PolyhedronNode* nodeNeighbour = (PolyhedronNode*)this->gcs.getNode(it->nodeKey);
			//Vertically concatenate the polyhedra constraints
			Eigen::MatrixXd bigA(convexSet->A.rows() + nodeNeighbour->polyhedron.A.rows(), convexSet->A.cols());
			Eigen::VectorXd bigB(convexSet->A.rows() + nodeNeighbour->polyhedron.A.rows());
			bigA << convexSet->A, nodeNeighbour->polyhedron.A;
			bigB << convexSet->b, nodeNeighbour->polyhedron.b;
			PolyhedronNode* edgeNode = new PolyhedronNode(bigA, bigB);
			//int nodeEGCSKey = this->egcs.addNode(edgeNode);
			int nodeEGCSKey = this->navGraph.addNode(edgeNode);
			neighbourNodeEGCSKeys.push_back(nodeEGCSKey);
			//Set the EGC2 to GCS map
			this->egcs2gcs.emplace(std::make_pair(nodeKey, it->nodeKey), nodeEGCSKey);
			//this->egcs2gcs.emplace(std::make_pair(it->nodeKey,nodeKey), nodeEGCSKey);
			//Set neighbours to EGCS
			for (std::map<std::pair<int, int>, int>::iterator itEdgeMap = this->egcs2gcs.begin(); itEdgeMap != this->egcs2gcs.end(); itEdgeMap++)
			{
				if ((itEdgeMap->first.first == it->nodeKey) || (itEdgeMap->first.second == it->nodeKey))
				{
					int foundEGCSKey = itEdgeMap->second;
					if (foundEGCSKey != nodeEGCSKey)
					{
						//this->egcs.addEdge(nodeEGCSKey, foundEGCSKey);
						//this->egcs.addEdge(foundEGCSKey, nodeEGCSKey);
						this->navGraph.addEdge(nodeEGCSKey, foundEGCSKey);
						this->navGraph.addEdge(foundEGCSKey, nodeEGCSKey);
						//We need to add also connections between foundEGCSKey and their terminal EGCS Keys
						for (std::vector<std::pair<int, int>>::iterator itTerminal = terminalEGCSNodeKeys.begin(); itTerminal != terminalEGCSNodeKeys.end(); itTerminal++)
						{
							if (itTerminal->first == it->nodeKey)
							{
								//this->egcs.addEdge(nodeEGCSKey, itTerminal->second); //One-way connection
								this->navGraph.addEdge(nodeEGCSKey, itTerminal->second); //One-way connection
							}
						}
					}
				}
				if ((itEdgeMap->first.first == nodeKey) || (itEdgeMap->first.second == nodeKey))
					//if ((itEdgeMap->first.first == nodeKey) )
				{
					int foundEGCSKey = itEdgeMap->second;
					if (foundEGCSKey != nodeEGCSKey)
					{
						//this->egcs.addEdge(nodeEGCSKey, foundEGCSKey);
						//this->egcs.addEdge(foundEGCSKey, nodeEGCSKey);
						this->navGraph.addEdge(nodeEGCSKey, foundEGCSKey);
						this->navGraph.addEdge(foundEGCSKey, nodeEGCSKey);
						//We need to add also connections between foundEGCSKey and their terminal EGCS Keys
						/*for (std::vector<std::pair<int, int>>::iterator itTerminal = terminalEGCSNodeKeys.begin(); itTerminal != terminalEGCSNodeKeys.end(); itTerminal++)
						{
							if (itTerminal->first == nodeKey)
								this->egcs->addEdge(nodeEGCSKey, itTerminal->second); //One-way connection
						}*/
					}
				}
			}
			//Add possible connections to the start node if possible
			bool isStartInsideNeighbour = nodeNeighbour->polyhedron.isInside(qStartNode->point.p);
			if (isStartInsideNeighbour)
			{
				this->navGraph.addEdge(this->qStartNodeNavGraphKey, nodeEGCSKey);
				/*for (std::vector<std::pair<int, int>>::iterator itTerminalEGCS = terminalEGCSNodeKeys.begin(); itTerminalEGCS != terminalEGCSNodeKeys.end(); itTerminalEGCS++)
				{
					if (itTerminalEGCS->first == nodeKey)
					{
						this->navGraph.addEdge(this->qStartNodeNavGraphKey, itTerminalEGCS->second);
					}
				}*/
			}

			//In addition, add one-way connections with recently created EGCS terminal nodes
			for (std::vector<std::pair<int, int>>::iterator itTerminal = terminalEGCSNodeKeys.begin(); itTerminal != terminalEGCSNodeKeys.end(); itTerminal++)
			{
				if (itTerminal->first == nodeKey)
				{
					//this->egcs.addEdge(nodeEGCSKey, itTerminal->second); //One-way connection
					this->navGraph.addEdge(nodeEGCSKey, itTerminal->second); //One-way connection
				}
			}
		}
		//Check if the original seed is inside the generated convex set. If not, repeat
		if (convexSet->isInside(q))
			break;
		ellipsoid = Ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	}
	//convexSet->print();
	//Now they have generated, we can treat them as obstacles for nexts calls
	/*for (std::vector<int>::iterator it = newNodeKeys.begin(); it != newNodeKeys.end(); it++)
	{
		Node* node = gcs->getNode(*it);
		PolyhedronNode* polyNode=(PolyhedronNode*)node->getNodeData();
		polyNode->useShrinked = false;
	}*/
}

void AStarIRISConic::buildNavGraph(const Eigen::VectorXd& qstart, const Eigen::VectorXd& qtarget)
{
	//this->navGraph = this->egcs;
	//int qstartKey=this->gcs.findConvexSet(qstart);
	//TODO: clear navGraph first
	this->qStartNode = new PointNode(qstart);
	this->qTargetNode = new PointNode(qtarget);
	this->qStartNodeNavGraphKey = this->navGraph.addNode(this->qStartNode);
	this->qTargetNodeNavGraphKey = this->navGraph.addNode(this->qTargetNode);
	/*for (std::vector<std::pair<int, int>>::iterator it = terminalEGCSNodeKeys.begin(); it != terminalEGCSNodeKeys.end(); it++)
	{
		if (it->first == qstartKey)
		{
			this->navGraph.addEdge(this->qstartNodeNavGraphKey, it->second);
		}
		this->navGraph.addEdge(it->second, this->qtargetNodeNavGraphKey);
	}
	for (std::map<std::pair<int, int>, int>::iterator itEdgeMap = this->egcs2gcs.begin(); itEdgeMap != this->egcs2gcs.end(); itEdgeMap++)
	{
		if ((itEdgeMap->first.first == qstartKey) || (itEdgeMap->first.second == qstartKey))
		{
			this->navGraph.addEdge(this->qstartNodeNavGraphKey, itEdgeMap->second);
		}
	}*/
}

void AStarIRISConic::terminalCosts(Eigen::VectorXd& hCosts, Eigen::MatrixXd& pClosest)
{
	Node* node = this->navGraph.getNode(this->qTargetNodeNavGraphKey);
	PointNode* targetNode = (PointNode*)node->getNodeData();
	hCosts.resize(terminalEGCSNodeKeys.size());
	pClosest.resize(targetNode->point.p.rows(),terminalEGCSNodeKeys.size());
	int i = 0;
	for (std::vector<std::pair<int, int>>::iterator it = terminalEGCSNodeKeys.begin(); it != terminalEGCSNodeKeys.end(); it++)
	{
		it->second;
		node=this->navGraph.getNode(it->second);
		PolyhedronTerminalNode* polyNode = (PolyhedronTerminalNode*)node->getNodeData();
		Eigen::VectorXd pOut;
		hCosts(i)=(polyNode->polyhedron.closestPoint(targetNode->point.p, pOut));
		pClosest.col(i) = pOut; 
		i++;
	}
}

void AStarIRISConic::expandTerminalNode(const int& terminalNodeKey)
{
	Node* node = this->navGraph.getNode(terminalNodeKey);
	PolyhedronTerminalNode* terminalNode = (PolyhedronTerminalNode*)node->getNodeData();
	if (terminalNode->expanded)
		return;
	int trials = 0;
	Eigen::VectorXd seed;
	while (trials<params.maxTrialsTerminal)
	{
		if (terminalNode->generateRandomSeed(seed,params.IRISParams.tol,params.maxTrialsTerminal))
		{
			
			if (!this->gcs.contains(seed, params.IRISParams.tol))  //This is not in Matlab code!!
			{
				if (this->_cObs->isFree(seed, params.IRISParams.tol))
				{
					trials = 0;
					//std::cout << "seed=" << std::endl;
					//std::cout << seed << std::endl;
					this->addConvexSets(seed);
				}
				else
				{
					trials++;
				}
			}
			else
			{
				trials++;
			}
		}
		else
		{
			trials++;
		}
	}
	terminalNode->expanded = true;
	//this->egcs.removeNode(terminalNodeKey);
	this->navGraph.removeNode(terminalNodeKey);
	for (std::vector<std::pair<int, int>>::iterator it = this->terminalEGCSNodeKeys.begin(); it != this->terminalEGCSNodeKeys.end(); it++)
	{
		if (it->second == terminalNodeKey)
		{
			this->terminalEGCSNodeKeys.erase(it);
			break;
		}
	}
}

void AStarIRISConic::doPhase1()
{
	this->addConvexSets(qStartNode->point.p);
	while (!this->gcs.contains(qTargetNode->point.p))
	{
		ConvexRelaxationMinDistance solver(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solver.setTask();
		solver.solve();
		//Select the terminal node connected to the target that generated the feasible solution (the size of nodeKeys is always at least 2, because it includes the start and target nodes)
		int terminalNodeKey=solver.feasibleSolution.nodeKeys[solver.feasibleSolution.nodeKeys.size()-2];
		//std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;
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
		std::cout << solver.relaxedSolution.cost << std::endl;
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
		std::cout << solver.feasibleSolution.cost << std::endl;*/
	}
}


AStarIRISParams_t AStarIRISConic::getDefaultAStarIRISParams()
{
	AStarIRISParams_t params = { IRISConic::getDefaultIRISParams(),99};
	return params;
}