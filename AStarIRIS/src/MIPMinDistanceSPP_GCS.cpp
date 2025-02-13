#include "MIPMinDistanceSPP_GCS.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PointNode.h"
#include "fusion.h"
#include <vector>
#include "EigenUtils.h"

MIPMinDistanceSPP_GCS::MIPMinDistanceSPP_GCS(Graph* g, const int& startKey, const int& targetKey) : ConvexRelaxationMinDistanceSPP_GCS(g,startKey, targetKey)
{
}

void MIPMinDistanceSPP_GCS::allocateVariables(Model::t& M)
{
	int nE = this->g->numEdges;
	this->l = M->variable("l", nE, Domain::greaterThan(0.0));
	this->z = M->variable("z", Set::make(this->n * this->N, nE), Domain::unbounded());
	this->y = M->variable("y", nE, Domain::integral(Domain::inRange(0, 1)));
}

/*void MinDistanceSolver::computeFeasibleSolution(const int& maxIters, const bool& simplifyGraph)
{
	int nextKey = startKey;

	while (nextKey != targetKey)
	{
		this->optimalPath.nodeKeys.push_back(nextKey);
		std::vector<int> outEdges = this->g->findOutEdges(nextKey);
		double max_y = 0.;
		int best_edge = -1;
		for (std::vector<int>::iterator it = outEdges.begin(); it != outEdges.end(); it++)
		{
			double current_y = this->MIPSolution.y(*it);
			if (current_y > max_y)
			{
				max_y = current_y;
				best_edge = *it;
			}
		}
		nextKey = this->g->getEdge(best_edge).second;
		this->optimalPath.edgeKeys.push_back(best_edge);
	}
	this->optimalPath.nodeKeys.push_back(nextKey);

	this->feasibleSolution.x.resize(this->optimalPath.nodeKeys.size(), this->MIPSolution.x.cols());
	int i = 0;
	this->feasibleSolution.cost = 0.;
	std::vector<int> nodeKeys = this->g->getNodeKeys();
	for (std::vector<int>::iterator it = this->optimalPath.nodeKeys.begin(); it != this->optimalPath.nodeKeys.end(); it++)
	{
		std::vector<int>::iterator it1 = std::find(nodeKeys.begin(), nodeKeys.end(), *it);
		int j = it1- nodeKeys.begin();
		this->feasibleSolution.x.row(i) = this->MIPSolution.x.row(j);
		i++;
	}
	i = 0;
	for (std::vector<int>::iterator it = this->optimalPath.edgeKeys.begin(); it != this->optimalPath.edgeKeys.end(); it++)
	{
		int j = *it;
		this->feasibleSolution.cost += this->MIPSolution.l(j);
		i++;
	}
}*/