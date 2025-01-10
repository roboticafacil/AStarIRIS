#include "ExpandableIRISConic.h"
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


ExpandableIRISConic::ExpandableIRISConic(CObsConic& cObs, const ExpandableIRISParams_t& params) : IRISConic(cObs, params.IRISParams), params(params)
{

}

ExpandableIRISConic::~ExpandableIRISConic()
{
	std::vector<int> keys=this->gcs.getNodeKeys();
	for (std::vector<int>::iterator it = keys.begin();it!=keys.end();it++)
	{
		Node* node=this->gcs.getNode(*it);
		this->gcs.removeNode(*it);
		delete(node);
	}
	keys = this->navGraph.getNodeKeys();
	for (std::vector<int>::iterator it = keys.begin(); it != keys.end(); it++)
	{
		Node* node = this->navGraph.getNode(*it);
		this->navGraph.removeNode(*it);
		delete(node);
	}
}

void ExpandableIRISConic::addConvexSets(const Eigen::VectorXd& q)
{
	Ellipsoid ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	Polyhedron* convexSet = new Polyhedron(q.rows());
	//std::vector<IRISNeighbour_t> neighbours;
	std::vector<int> neighbourKeys;
	int count = 0;
	int countTerminal = 0;
	int countIntersection = 0;
	std::vector<int> newNodeKeys;
	std::vector<int> newTerminalNodeKeys;
	std::vector<int> newIntersectionNodeKeys;
	//convexSet.allocateClosestPointEllipsoidSolver();
	//convexSet.allocateInscribedEllipsoidSolver();
	while (true)
	{
		//IRISConic::computeConvexSet(ellipsoid, *convexSet, neighbours);
		IRISConic::computeConvexSet(ellipsoid, *convexSet,neighbourKeys);
		//std::cout << "Convex set added to GCS " << std::endl;
		
		//Add the convexSet to GCS
		Eigen::VectorXd bShrinked = (1. - this->params.IRISParams.shrinkFactor) * convexSet->b + this->params.IRISParams.shrinkFactor * (convexSet->A * ellipsoid.getCentroid());
		PolyhedronNode* node = new PolyhedronNode(convexSet->A, convexSet->b, bShrinked);
		bool isStartInside = convexSet->isInside(qStartNode->point.p);
		bool isTargetInside = convexSet->isInside(qTargetNode->point.p);
		int nodeKey = this->gcs.addNode(node);
		//if (nodeKey == 21)
		//	std::cout << "The next one seems problematic" << std::endl;
		newNodeKeys.push_back(nodeKey);
		//std::cout << "Convex Set " << nodeKey << std::endl;
		//convexSet->print();
		//std::cout << "Ellipsoid" << nodeKey << std::endl;
		//ellipsoid.print();
		count++;
		//newNodeKeys.push_back(nodeKey);
		//Add terminal nodes to EGCS
		//std::cout << "Ellipsoid" << std::endl;
		//ellipsoid.print();
		for (int i = 0; i < convexSet->A.rows(); i++)
		{
			Eigen::VectorXd ai = -convexSet->A.row(i);
			double bi = -(convexSet->b(i) - params.IRISParams.tol);
			//std::cout << "ai=" << std::endl;
			//std::cout << ai << std::endl;
			//std::cout << "bi=" << bi << std::endl;
			//if (!ellipsoid.isInsideSeparatingHyperplane(ai,bi,0.))
			//	continue;
			Eigen::MatrixXd facetA(convexSet->A.rows() + 1, convexSet->A.cols());
			Eigen::VectorXd facetB(convexSet->A.rows() + 1);
			//Eigen::VectorXd facetB1(convexSet->A.rows() + 1);
			Eigen::VectorXd b = convexSet->b;
			b(i) = b(i) + params.IRISParams.tol; //Expand slightly the original polyhedron in the facet direction
			facetA << convexSet->A, -convexSet->A.row(i);
			facetB << b, -convexSet->b(i);
			//facetB1 << convexSet->b, -convexSet->b(i) + params.IRISParams.tol;
			//std::cout << "Facet A" << std::endl;
			//std::cout << facetA << std::endl;
			//std::cout << "facet B" << std::endl;
			//std::cout << facetB << std::endl;
			//std::cout << "facet B1" << std::endl;
			//std::cout << facetB1 << std::endl;
			//int terminalEGCSKey = this->egcs.addNode(new PolyhedronTerminalNode(facetA, facetB,i));
			int terminalNavGraphKey = this->navGraph.addNode(new PolyhedronTerminalNode(facetA, facetB, i));
			newTerminalNodeKeys.push_back(terminalNavGraphKey);
			countTerminal++;
			terminalNodeKeys.push_back(std::make_pair(nodeKey, terminalNavGraphKey));
			if (isStartInside)
				this->navGraph.addEdge(this->qStartNodeNavGraphKey, terminalNavGraphKey);
			this->navGraph.addEdge(terminalNavGraphKey, this->qTargetNodeNavGraphKey);
		}
		//Add neighbours to GCS, create EGCS node, set the EGCS to GCS map and set neighbours of EGCS node
		std::vector<int> neighbourNodeEGCSKeys;
		for (std::vector<int>::iterator it = neighbourKeys.begin(); it != neighbourKeys.end(); it++)
		{
			//Add neighbours to GCS
			int neighbourKey = *it;
			this->gcs.addEdge(nodeKey, neighbourKey);
			//this->gcs->addEdge(it->nodeKey,nodeKey);  //Maybe this is not necessary
			//Create EGCS node
			PolyhedronNode* nodeNeighbour = (PolyhedronNode*)this->gcs.getNode(neighbourKey);
			//Vertically concatenate the polyhedra constraints
			Eigen::MatrixXd bigA(convexSet->A.rows() + nodeNeighbour->polyhedron.A.rows(), convexSet->A.cols());
			Eigen::VectorXd bigB(convexSet->A.rows() + nodeNeighbour->polyhedron.A.rows());
			bigA << convexSet->A, nodeNeighbour->polyhedron.A;
			bigB << convexSet->b, nodeNeighbour->polyhedron.b;
			PolyhedronNode* edgeNode = new PolyhedronNode(bigA, bigB);
			//int nodeEGCSKey = this->egcs.addNode(edgeNode);
			int nodeIntersectionKey = this->navGraph.addNode(edgeNode);
			newIntersectionNodeKeys.push_back(nodeIntersectionKey);
			countIntersection++;
			neighbourNodeEGCSKeys.push_back(nodeIntersectionKey);
			//Set the EGC2 to GCS map
			this->navGraph2gcs.emplace(std::make_pair(nodeKey, neighbourKey), nodeIntersectionKey);
			if (isTargetInside)
				this->navGraph.addEdge(nodeIntersectionKey, this->qTargetNodeNavGraphKey);
			bool isInsideNeighbour=nodeNeighbour->polyhedron.isInside(qTargetNode->point.p);
			if (isInsideNeighbour)
				this->navGraph.addEdge(nodeIntersectionKey, this->qTargetNodeNavGraphKey);
			//this->egcs2gcs.emplace(std::make_pair(it->nodeKey,nodeKey), nodeEGCSKey);
			//Set neighbours to EGCS
			for (std::map<std::pair<int, int>, int>::iterator itEdgeMap = this->navGraph2gcs.begin(); itEdgeMap != this->navGraph2gcs.end(); itEdgeMap++)
			{
				if ((itEdgeMap->first.first == neighbourKey) || (itEdgeMap->first.second == neighbourKey))
				{
					int foundEGCSKey = itEdgeMap->second;
					if (foundEGCSKey != nodeIntersectionKey)
					{
						//this->egcs.addEdge(nodeEGCSKey, foundEGCSKey);
						//this->egcs.addEdge(foundEGCSKey, nodeEGCSKey);
						this->navGraph.addEdge(nodeIntersectionKey, foundEGCSKey);
						this->navGraph.addEdge(foundEGCSKey, nodeIntersectionKey);
						//We need to add also connections between foundEGCSKey and their terminal EGCS Keys
						for (std::vector<std::pair<int, int>>::iterator itTerminal = terminalNodeKeys.begin(); itTerminal != terminalNodeKeys.end(); itTerminal++)
						{
							if (itTerminal->first == neighbourKey)
							{
								this->navGraph.addEdge(nodeIntersectionKey, itTerminal->second); //One-way connection
							}
						}
					}
				}
				if ((itEdgeMap->first.first == nodeKey) || (itEdgeMap->first.second == nodeKey))
					//if ((itEdgeMap->first.first == nodeKey) )
				{
					int foundEGCSKey = itEdgeMap->second;
					if (foundEGCSKey != nodeIntersectionKey)
					{
						//this->egcs.addEdge(nodeEGCSKey, foundEGCSKey);
						//this->egcs.addEdge(foundEGCSKey, nodeEGCSKey);
						this->navGraph.addEdge(nodeIntersectionKey, foundEGCSKey);
						this->navGraph.addEdge(foundEGCSKey, nodeIntersectionKey);
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
				this->navGraph.addEdge(this->qStartNodeNavGraphKey, nodeIntersectionKey);
				/*for (std::vector<std::pair<int, int>>::iterator itTerminalEGCS = terminalEGCSNodeKeys.begin(); itTerminalEGCS != terminalEGCSNodeKeys.end(); itTerminalEGCS++)
				{
					if (itTerminalEGCS->first == nodeKey)
					{
						this->navGraph.addEdge(this->qStartNodeNavGraphKey, itTerminalEGCS->second);
					}
				}*/
			}

			//In addition, add one-way connections with recently created EGCS terminal nodes
			for (std::vector<std::pair<int, int>>::iterator itTerminal = terminalNodeKeys.begin(); itTerminal != terminalNodeKeys.end(); itTerminal++)
			{
				if (itTerminal->first == nodeKey)
				{
					this->navGraph.addEdge(nodeIntersectionKey, itTerminal->second); //One-way connection
				}
			}
		}
		//Check if the original seed is inside the generated convex set. If not, repeat
		if (convexSet->isInside(q))
			break;
		ellipsoid = Ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	}
	if (count > 0)
	{
		std::cout << count << " new convex sets added to GCS: ";
		for (std::vector<int>::iterator it = newNodeKeys.begin(); it != newNodeKeys.end(); it++)
			std::cout << *it << " ";
		std::cout << std::endl;
		std::cout << countIntersection << " new Intersection convex sets added to navGraph: ";
		for (std::vector<int>::iterator it = newIntersectionNodeKeys.begin(); it != newIntersectionNodeKeys.end(); it++)
			std::cout << *it << " ";
		std::cout << std::endl;
		std::cout << countTerminal << " new Terminal convex sets added to navGraph: ";
		for (std::vector<int>::iterator it = newTerminalNodeKeys.begin(); it != newTerminalNodeKeys.end(); it++)
			std::cout << *it << " ";
		std::cout << std::endl; 
	}
	//convexSet->print();
	//Now they have generated, we can treat them as obstacles for nexts calls
	/*for (std::vector<int>::iterator it = newNodeKeys.begin(); it != newNodeKeys.end(); it++)
	{
		Node* node = gcs->getNode(*it);
		PolyhedronNode* polyNode=(PolyhedronNode*)node->getNodeData();
		polyNode->useShrinked = false;
	}*/
	delete(convexSet);
}

void ExpandableIRISConic::buildNavGraph(const Eigen::VectorXd& qstart, const Eigen::VectorXd& qtarget)
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

void ExpandableIRISConic::terminalCosts(Eigen::VectorXd& hCosts, Eigen::MatrixXd& pClosest)
{
	Node* node = this->navGraph.getNode(this->qTargetNodeNavGraphKey);
	PointNode* targetNode = (PointNode*)node->getNodeData();
	hCosts.resize(terminalNodeKeys.size());
	pClosest.resize(targetNode->point.p.rows(), terminalNodeKeys.size());
	int i = 0;
	for (std::vector<std::pair<int, int>>::iterator it = terminalNodeKeys.begin(); it != terminalNodeKeys.end(); it++)
	{
		it->second;
		node = this->navGraph.getNode(it->second);
		PolyhedronTerminalNode* polyNode = (PolyhedronTerminalNode*)node->getNodeData();
		Eigen::VectorXd pOut;
		hCosts(i) = (polyNode->polyhedron.closestPoint(targetNode->point.p, pOut));
		pClosest.col(i) = pOut;
		i++;
	}
}

void ExpandableIRISConic::expandTerminalNode(const int& terminalNodeKey)
{
	Node* node = this->navGraph.getNode(terminalNodeKey);
	if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
	{
		PolyhedronTerminalNode* terminalNode = (PolyhedronTerminalNode*)node->getNodeData();
		if (terminalNode->expanded)
		{
			delete(terminalNode);
			return;
		}
		int trials = 0;
		Eigen::VectorXd seed;
		while (trials < params.maxTrialsTerminal)
		{
			if (terminalNode->generateRandomSeed(seed, params.IRISParams.tol, params.maxTrialsTerminal))
			{
				//std::cout << "seed=" << std::endl;
				//std::cout << seed << std::endl;
				//if (!this->gcs.contains(seed, params.IRISParams.tol))  //This is not in Matlab code!!
				if (!this->gcs.contains(seed, 0.))  //This is not in Matlab code!!
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
		//this->egcs.removeNode(terminalNodeKey);
		/*this->navGraph.removeNode(terminalNodeKey);
		for (std::vector<std::pair<int, int>>::iterator it = this->terminalEGCSNodeKeys.begin(); it != this->terminalEGCSNodeKeys.end(); it++)
		{
			if (it->second == terminalNodeKey)
			{
				this->terminalEGCSNodeKeys.erase(it);
				break;
			}
		}
		delete(node);
		terminalNode->expanded = true;*/
		this->removeNode(terminalNodeKey);
	}
	else
		std::cout << "Warning:: Trying to expand a node that is NOT a terminal node" << std::endl;
}

NavGraph ExpandableIRISConic::getGraphWithoutTerminalConnections()
{
	NavGraph simplifiedNavGraph = this->navGraph;
	std::vector<int> nodeKeys = simplifiedNavGraph.getNodeKeys();
	std::vector<Edge> edges = simplifiedNavGraph.getEdges();
	std::vector<Edge> deletedEdges;
	for (std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); )
	{
		if (it->second == qTargetNodeNavGraphKey)
		{
			Node* node = (Node*)this->navGraph.getNode(it->first);
			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				//deletedEdges.push_back(*it);
				//edges.erase(it);
				simplifiedNavGraph.removeNode(it->first);
				it++;
			}
			else
				it++;
		}
		else
			it++;
	}
	//std::cout << "Deleted Edges " << std::endl;
	//for (std::vector<Edge>::iterator it = deletedEdges.begin(); it != deletedEdges.end(); it++)
	//	std::cout << "Deleted edge: " << it->first << "->" << it->second << std::endl;
	//simplifiedNavGraph.setEdges(edges);
	return simplifiedNavGraph;
}

void ExpandableIRISConic::removeNode(const int key)
{
	Node* node = this->navGraph.getNode(key);
	this->navGraph.removeNode(key);
	for (std::vector<std::pair<int, int>>::iterator it = this->terminalNodeKeys.begin(); it != this->terminalNodeKeys.end(); it++)
	{
		if (it->second == key)
		{
			this->terminalNodeKeys.erase(it);
			break;
		}
	}
	delete(node);
}


ExpandableIRISParams_t ExpandableIRISConic::getDefaultExpandableIRISParams()
{
	ExpandableIRISParams_t params = { IRISConic::getDefaultIRISParams(),99, 0 };
	return params;
}