#include "Solver.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "fusion.h"
#include <vector>
#include "EigenNdArray.h"

Solver::Solver(const Graph& g, const int& startKey, const int& targetKey) : g(g), nN(g.numNodes), nE(g.numEdges), startKey(startKey), targetKey(targetKey)
{
	PointNode* qs = (PointNode*)this->g.getNode(startKey);
	Eigen::MatrixXd qsEigen(qs->point.p);
	this->n = qs->point.p.rows();
	qstart = Eigen2NdArray(qsEigen);
	PointNode* qt = (PointNode*)this->g.getNode(targetKey);
	Eigen::MatrixXd qtEigen(qt->point.p);
	qtarget = Eigen2NdArray(qtEigen);
}

Solver::~Solver()
{
	if (solverAllocated)
	{
		solverAllocated = false;
		this->M->dispose();
	}
}

bool Solver::optimalPathContainsTerminalNodes()
{
	for (std::vector<int>::iterator it = this->optimalPath.nodeKeys.begin(); it != this->optimalPath.nodeKeys.end(); it++)
	{
		Node* node = this->g.getNode(*it);
		if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
		{
			return true;
		}
	}
	return false;
}

void Solver::setGraph(Graph& g)
{
	this->g = g;
}