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

int ExpandableIRISConic::addConvexSet(const Eigen::VectorXd& q)
{
	Ellipsoid ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	Polyhedron convexSet = Polyhedron(q.rows());
	std::vector<int> neighbourKeys;
	int count = 0;
	int countTerminal = 0;
	int countIntersection = 0;
	std::vector<int> newTerminalNodeKeys;
	std::vector<int> newIntersectionNodeKeys;

	this->computeConvexSet(ellipsoid, convexSet, neighbourKeys);
	//Add the convexSet to GCS
	//ellipsoid.print();
	//convexSet.print();
	//std::cout << "q: " << std::endl;
	//std::cout << q << std::endl;
	Eigen::VectorXd bShrinked = (1. - this->params.IRISParams.shrinkFactor) * convexSet.b + this->params.IRISParams.shrinkFactor * (convexSet.A * ellipsoid.getCentroid());
	PolyhedronNode* node = new PolyhedronNode(convexSet.A, convexSet.b, bShrinked);
	bool isStartInside = convexSet.isInside(qStartNode->point.p);
	bool isTargetInside = convexSet.isInside(qTargetNode->point.p);
	int nodeKey = this->gcs.addNode(node);
	//std::cout << "Convex Set " << nodeKey << std::endl;
	//convexSet.print();
	//std::cout << "bShrinked: " << std::endl;
	//std::cout << bShrinked << std::endl;
	//std::cout << "Ellipsoid" << nodeKey << std::endl;
	//ellipsoid.print();
	count++;
	//Add terminal nodes to EGCS
	//std::cout << "Ellipsoid" << std::endl;
	//ellipsoid.print();
	for (int i = 0; i < convexSet.A.rows(); i++)
	{
		Eigen::VectorXd ai = -convexSet.A.row(i);
		double bi = -(convexSet.b(i) - params.IRISParams.tol);
		//std::cout << "ai=" << std::endl;
		//std::cout << ai << std::endl;
		//std::cout << "bi=" << bi << std::endl;
		//if (!ellipsoid.isInsideSeparatingHyperplane(ai,bi,0.))
		//	continue;
		Eigen::MatrixXd facetA(convexSet.A.rows() + 1, convexSet.A.cols());
		Eigen::VectorXd facetB(convexSet.A.rows() + 1);
		//Eigen::VectorXd facetB1(convexSet.A.rows() + 1);
		Eigen::VectorXd b = convexSet.b;
		b(i) = b(i) + params.IRISParams.tol; //Expand slightly the original polyhedron in the facet direction
		facetA << convexSet.A, -convexSet.A.row(i);
		facetB << b, -convexSet.b(i);
		//facetB1 << convexSet.b, -convexSet.b(i) + params.IRISParams.tol;
		//std::cout << "Facet A" << std::endl;
		//std::cout << facetA << std::endl;
		//std::cout << "facet B" << std::endl;
		//std::cout << facetB << std::endl;
		//std::cout << "facet B1" << std::endl;
		//std::cout << facetB1 << std::endl;
		int terminalNavGraphKey = this->navGraph.addNode(new PolyhedronTerminalNode(facetA, facetB, i));
		newTerminalNodeKeys.push_back(terminalNavGraphKey);
		countTerminal++;
		terminalNodeKeys.push_back(std::make_pair(nodeKey, terminalNavGraphKey));
		if (isStartInside)
			this->navGraph.addEdge(this->qStartNodeNavGraphKey, terminalNavGraphKey);
		this->navGraph.addEdge(terminalNavGraphKey, this->qTargetNodeNavGraphKey);
	}
	//Add neighbours to GCS, create navGraph node, set the navGraph to GCS map and set neighbours of navGraph node
	std::vector<int> neighbourNodeNavGraphKeys;
	for (std::vector<int>::iterator it = neighbourKeys.begin(); it != neighbourKeys.end(); it++)
	{
		//Add neighbours to GCS
		int neighbourKey = *it;
		this->gcs.addEdge(nodeKey, neighbourKey);
		//this->gcs->addEdge(it->nodeKey,nodeKey);  //Maybe this is not necessary
		//Create navGraph node
		PolyhedronNode* nodeNeighbour = (PolyhedronNode*)this->gcs.getNode(neighbourKey);
		//Vertically concatenate the polyhedra constraints
		Eigen::MatrixXd bigA(convexSet.A.rows() + nodeNeighbour->polyhedron.A.rows(), convexSet.A.cols());
		Eigen::VectorXd bigB(convexSet.A.rows() + nodeNeighbour->polyhedron.A.rows());
		bigA << convexSet.A, nodeNeighbour->polyhedron.A;
		bigB << convexSet.b, nodeNeighbour->polyhedron.b;
		PolyhedronNode* edgeNode = new PolyhedronNode(bigA, bigB);
		int nodeIntersectionKey = this->navGraph.addNode(edgeNode);
		newIntersectionNodeKeys.push_back(nodeIntersectionKey);
		countIntersection++;
		neighbourNodeNavGraphKeys.push_back(nodeIntersectionKey);
		//Set the navGraph to GCS map
		this->navGraph2gcs.emplace(std::make_pair(nodeKey, neighbourKey), nodeIntersectionKey);
		if (isTargetInside)
			this->navGraph.addEdge(nodeIntersectionKey, this->qTargetNodeNavGraphKey);
		bool isInsideNeighbour = nodeNeighbour->polyhedron.isInside(qTargetNode->point.p);
		if (isInsideNeighbour)
			this->navGraph.addEdge(nodeIntersectionKey, this->qTargetNodeNavGraphKey);
		//Set neighbours to navGraph node
		for (std::map<std::pair<int, int>, int>::iterator itEdgeMap = this->navGraph2gcs.begin(); itEdgeMap != this->navGraph2gcs.end(); itEdgeMap++)
		{
			if ((itEdgeMap->first.first == neighbourKey) || (itEdgeMap->first.second == neighbourKey))
			{
				int foundNavGraphKey = itEdgeMap->second;
				if (foundNavGraphKey != nodeIntersectionKey)
				{
					this->navGraph.addEdge(nodeIntersectionKey, foundNavGraphKey);
					this->navGraph.addEdge(foundNavGraphKey, nodeIntersectionKey);
					//We need to add also connections between foundNavGraphKey and their terminal navGraph Keys
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
			{
				int foundNavGraphKey = itEdgeMap->second;
				if (foundNavGraphKey != nodeIntersectionKey)
				{
					this->navGraph.addEdge(nodeIntersectionKey, foundNavGraphKey);
					this->navGraph.addEdge(foundNavGraphKey, nodeIntersectionKey);
				}
			}
		}
		//Add possible connections to the start node if possible
		bool isStartInsideNeighbour = nodeNeighbour->polyhedron.isInside(qStartNode->point.p);
		if (isStartInsideNeighbour)
		{
			this->navGraph.addEdge(this->qStartNodeNavGraphKey, nodeIntersectionKey);
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
	//if (convexSet->isInside(q))
	//	break;
	ellipsoid = Ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);

	std::cout << "New convex set added to GCS: " << nodeKey << std::endl;
	std::cout << countIntersection << " new Intersection convex sets added to navGraph: ";
	for (std::vector<int>::iterator it = newIntersectionNodeKeys.begin(); it != newIntersectionNodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << countTerminal << " new Terminal convex sets added to navGraph: ";
	for (std::vector<int>::iterator it = newTerminalNodeKeys.begin(); it != newTerminalNodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	return nodeKey;
}

std::vector<int> ExpandableIRISConic::addConvexSets(const Eigen::VectorXd& q)
{
	std::vector<int> convexSetNodeKeys;
	while (true)
	{
		int convexSetNodeKey = this->addConvexSet(q);
		convexSetNodeKeys.push_back(convexSetNodeKey);
		PolyhedronNode* node=(PolyhedronNode*)this->gcs.getNode(convexSetNodeKey);
		//Check if the original seed is inside the generated convex set. If not, repeat
		if (node->polyhedron.isInside(q))
			break;
	}
	return convexSetNodeKeys;
}

void ExpandableIRISConic::computeConvexSet(Ellipsoid& ellipsoid, Polyhedron& convexSet, std::vector<int>& neighbourKeys)
{
	Ellipsoid newEllipsoid(ellipsoid);
	double detC = 0.;
	Eigen::VectorXd q0 = ellipsoid.getCentroid();
	Polyhedron newConvexSet(convexSet.n);
	while (true)
	{
		//std::cout << "Computing separating hyperplanes with ellipsoid " << std::endl;
		//newEllipsoid.print();
		separatingHyperplanes(newEllipsoid, newConvexSet);
		//std::cout << "IRIS: new convex set computed" << std::endl;
		//newConvexSet.print();
		if (!newConvexSet.isInside(q0))
		{
			//std::cout << "IRIS: seed " << std::endl;
			//std::cout << q0 << std::endl;
			break;
		}
		newEllipsoid = newConvexSet.inscribedEllipsoid();
		//newEllipsoid.print();
		double detNewC = newEllipsoid.C.determinant();
		if (((detNewC - detC) / detC) < params.IRISParams.tolConvexSetConvergence)
		{
			if ((params.IRISParams.seperatingHyperplaneAligned) && (detNewC < detC))
			{
				break;
			}
			convexSet = newConvexSet;
			ellipsoid = newEllipsoid;
			break;
		}
		convexSet = newConvexSet;
		ellipsoid = newEllipsoid;
		detC = detNewC;
	}
	convexSet.removeConstraints();
	neighbourKeys.clear();
	std::vector<int> nodeKeys = this->gcs.getNodeKeys();
	for (std::vector<int>::iterator itNode = nodeKeys.begin(); itNode != nodeKeys.end(); itNode++)
	{
		Node* existingNode = this->gcs.getNode(*itNode);
		PolyhedronNode* existingPolyNode = (PolyhedronNode*)existingNode;
		if (Polyhedron::intersect(convexSet, existingPolyNode->polyhedron))
			neighbourKeys.push_back(*itNode);
	}
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
			//std::cout << "Terminal node polyhedron" << std::endl;
			//terminalNode->polyhedron.print();
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
						this->addConvexSet(seed);
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