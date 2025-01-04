#include "ConvexRelaxation.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PointNode.h"
#include "fusion.h"
#include <vector>
#include "EigenNdArray.h"

ConvexRelaxation::ConvexRelaxation(const Graph& g, const int& startKey, const int& targetKey) : g(g), nN(g.numNodes), nE(g.numEdges), startKey(startKey), targetKey(targetKey)
{
	PointNode* qs = (PointNode*)this->g.getNode(startKey);
	Eigen::MatrixXd qsEigen(qs->point.p);
	this->n = qs->point.p.rows();
	qstart = Eigen2NdArray(qsEigen);
	PointNode* qt = (PointNode*)this->g.getNode(targetKey);
	Eigen::MatrixXd qtEigen(qt->point.p);
	qtarget = Eigen2NdArray(qtEigen);
}

void ConvexRelaxation::computeFeasibleSolution()
{	
	int nextKey = startKey;
	
	while (nextKey != targetKey)
	{
		this->feasibleSolution.nodeKeys.push_back(nextKey);
		std::vector<int> outEdges = this->g.findOutEdges(nextKey);
		double max_y = 0.;
		int best_edge = -1;
		for (std::vector<int>::iterator it = outEdges.begin(); it != outEdges.end(); it++)
		{
			double current_y = this->relaxedSolution.y(*it);
			
			if (current_y > max_y)
			{
				int current_key = this->g.getEdge(*it).second;
				bool not_found = true;
				for (std::vector<int>::iterator itNodeKey = this->feasibleSolution.nodeKeys.begin(); itNodeKey != this->feasibleSolution.nodeKeys.end(); itNodeKey++)
				{
					if (*itNodeKey == current_key)
					{
						not_found = false;
						break;
					}
				}
				if (not_found)
				{
					max_y = current_y;
					best_edge = *it;
				}
			}
		}
		nextKey=this->g.getEdge(best_edge).second;
		this->feasibleSolution.edgeKeys.push_back(best_edge);
	}
	this->feasibleSolution.nodeKeys.push_back(nextKey);
	
	this->feasibleSolution.x.resize(this->feasibleSolution.edgeKeys.size() + 1, this->relaxedSolution.x.cols());
	int i = 0;
	int j = 0;
	this->feasibleSolution.cost = 0.;
	for (std::vector<int>::iterator it = this->feasibleSolution.edgeKeys.begin(); it != this->feasibleSolution.edgeKeys.end(); it++)
	{
		j = *it;
		this->feasibleSolution.x.row(i)= this->relaxedSolution.z.row(j);
		this->feasibleSolution.cost += this->relaxedSolution.l(j);
		i++;
	}
	this->feasibleSolution.x.row(i) = this->relaxedSolution.p.row(j);
}

void ConvexRelaxation::solve()
{
	this->M->solve();
	this->relaxedSolution.status = this->M->getPrimalSolutionStatus();
	if (relaxedSolution.status == SolutionStatus::Optimal) {


		Eigen::VectorXd y=Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->M->getVariable("y")->level()->begin(), this->nE));
		Eigen::VectorXd l=Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->M->getVariable("l")->level()->begin(), this->nE));
		Eigen::VectorXd zvec=Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->M->getVariable("z")->level()->begin(), this->n * this->nE));
		Eigen::VectorXd pvec=Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->M->getVariable("p")->level()->begin(), this->n * this->nE));
		this->relaxedSolution.y = y;
		this->relaxedSolution.l = l;
		this->relaxedSolution.z = zvec.reshaped(this->nE, this->n).eval();
		this->relaxedSolution.p = pvec.reshaped(this->nE, this->n).eval();
		this->relaxedSolution.x = Eigen::MatrixXd(this->nN, this->n);
		std::vector<int> nodeKeys=this->g.getNodeKeys();
		for (int i = 0; i < this->nN; i++)
		{

			std::vector<int> inEdges = this->g.findInEdges(nodeKeys[i]);
			if (inEdges.size() > 0)
			{
				for (std::vector<int>::iterator it = inEdges.begin(); it != inEdges.end(); it++)
				{
					this->relaxedSolution.x.row(i) += this->relaxedSolution.y(*it) * this->relaxedSolution.p.row(*it);
				}
			}
			else
			{
				std::vector<int> outEdges = this->g.findOutEdges(nodeKeys[i]);
				if (outEdges.size() > 0)
				{
					for (std::vector<int>::iterator it = outEdges.begin(); it != outEdges.end(); it++)
					{
						this->relaxedSolution.x.row(i) += this->relaxedSolution.y(*it) * this->relaxedSolution.z.row(*it);
					}
				}
				else
				{
					std::cout << "Node " << i << " has no edges" << std::endl;
				}
			}
		}
		this->relaxedSolution.cost = this->M->primalObjValue();
	}
	else
	{
		std::cout << "Another solution status" << std::endl;
		std::cout << "Solution status: " << this->M->getPrimalSolutionStatus() << std::endl;
	}
	this->computeFeasibleSolution();
}