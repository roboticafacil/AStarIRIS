#include "AStarGCS.h"
#include "PolyhedronNode.h"
#include "PointNode.h"
#include "EllipsoidNode.h"
#include "PolyhedronObstacleCircularRobotNode.h"

AStarGCS::AStarGCS(const Graph& g, const std::vector<int>& nodeKeys, const int& startKey, const int& targetKey): graph(g), nodeKeys(nodeKeys), h(g.numNodes), nN(g.numNodes), startKey(startKey), targetKey(targetKey)
{
	PointNode* node = (PointNode*)graph.getNode(startKey);
	this->n = node->point.p.size();
	this->pClosest=Eigen::MatrixXd(n,nN);
	this->computeHeuristics();
}

void AStarGCS::computeHeuristics()
{
	PointNode* pTarget = (PointNode*)graph.getNode(targetKey);
	for (int i = 0; i < nN; i++)
	{
		Node* node = graph.getNode(this->nodeKeys[i]);
		ConicSet* conic = dynamic_cast<ConicSet*>((ConicSet*)node->getNodeData());
		Eigen::VectorXd pOut;
		this->h(i) = conic->closestPoint(pTarget->point.p, pOut);
		this->pClosest.col(i) = pOut;

		/*if (dynamic_cast<PolyhedronNode*>(node))
		{
			PolyhedronNode* polyNode = (PolyhedronNode*)node;
			Eigen::VectorXd pOut;
			this->h(i) = polyNode->polyhedron.closestPoint(pTarget->point.p,pOut);
			this->pClosest.col(i) = pOut;
		}
		else if (dynamic_cast<PointNode*>(node)) {
			PointNode* pNode = (PointNode*)node;
			this->h(i)=pNode->point.distance(pTarget->point.p);
			this->pClosest.col(i) = pNode->point.p;
		}
		else if (dynamic_cast<PolyhedronObstacleCircularRobotNode*>(node)) {
			PolyhedronObstacleCircularRobotNode* polyNode = (PolyhedronObstacleCircularRobotNode*)node;
			Eigen::VectorXd pOut;
			this->h(i) = polyNode->polyhedronObstalce.closestPoint(pTarget->point.p, pOut);
			this->pClosest.col(i) = pOut;
		}
		else if (dynamic_cast<EllipsoidNode*>(node)) {
			EllipsoidNode* eNode = (EllipsoidNode*)node;
			Eigen::VectorXd pOut;
			this->h(i) = eNode->ellipsoid.closestPoint(pTarget->point.p, pOut);
			this->pClosest.col(i) = pOut;
		}*/
	}
}

Eigen::VectorXd AStarGCS::getHeuristicCosts()
{
	return this->h;
}
