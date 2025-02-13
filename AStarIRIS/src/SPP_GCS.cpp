#include "SPP_GCS.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "fusion.h"
#include <vector>
#include "EigenUtils.h"

SPP_GCS::SPP_GCS(Graph* g, const int &N, const int& startKey, const int& targetKey) : g(g), N(N), startKey(startKey), targetKey(targetKey)
{
	PointNode* qs = (PointNode*)this->g->getNode(startKey);
	Eigen::MatrixXd qsEigen(qs->point.p);
	this->n = qs->point.p.rows();
	qstart = Eigen2NdArray(qsEigen);
	PointNode* qt = (PointNode*)this->g->getNode(targetKey);
	Eigen::MatrixXd qtEigen(qt->point.p);
	qtarget = Eigen2NdArray(qtEigen);
	/*int k = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < nE; j++)
			idx(i, j) = k++;
	}*/
}

bool SPP_GCS::optimalPathContainsTerminalNodes()
{
	for (std::vector<int>::iterator it = this->optimalPath.nodeKeys.begin(); it != this->optimalPath.nodeKeys.end(); it++)
	{
		Node* node = this->g->getNode(*it);
		if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
		{
			return true;
		}
	}
	return false;
}

void SPP_GCS::setGraph(Graph* g)
{
	this->g = g;
}