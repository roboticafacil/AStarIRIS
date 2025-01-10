#include "ConvexRelaxationSolver.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PointNode.h"
#include "fusion.h"
#include <vector>
#include "EigenNdArray.h"
#include <random>
#include <limits>

ConvexRelaxationSolver::ConvexRelaxationSolver(const Graph& g, const int& startKey, const int& targetKey) : Solver(g,startKey,targetKey)
{
}

bool ConvexRelaxationSolver::getMCPath(Path_t & path)
{
	int nextKey = startKey;
	//Path_t MCPath;
	std::random_device rd;

	path.nodeKeys.clear();
	path.edgeKeys.clear();
	
	int iters = 0;
	while (nextKey != targetKey)
	{
		path.nodeKeys.push_back(nextKey);
		std::vector<int> outEdges = this->g.findOutEdges(nextKey);
		std::vector<double> weights;
		for (std::vector<int>::iterator it = outEdges.begin(); it != outEdges.end(); it++)
		{
			weights.push_back(this->relaxedSolution.y(*it));
		}
		std::discrete_distribution<int> dist(weights.begin(), weights.end());

		
		int selected_edge;
		//bool found = true;
		
		//while (found)
		//{
			std::mt19937 rng(rd());
			selected_edge = outEdges[dist(rng)];
			nextKey = this->g.getEdge(selected_edge).second;
			std::vector<int>::iterator it=std::find(path.nodeKeys.begin(), path.nodeKeys.end(), nextKey);
			if (it < path.nodeKeys.end())
				return false;

			//found = false;
			/*for (std::vector<int>::iterator itNodeKey = path.nodeKeys.begin(); itNodeKey != path.nodeKeys.end(); itNodeKey++)
			{
				if (*itNodeKey == nextKey)
				{

					int nodesToDelete = path.nodeKeys.end() - (itNodeKey + 1);
					//std::cout << nodesToDelete << std::endl;
					path.nodeKeys.erase(itNodeKey + 1, path.nodeKeys.end());
					path.edgeKeys.erase(path.edgeKeys.end() - nodesToDelete, path.edgeKeys.end());

					found = true;
					


					//This is not working because redirecting it to another path might not reach the target
					//while (*itNodeKey == nextKey)
					//{
					//	selected_edge = dist(rng);
					//	nextKey = this->g.getEdge(selected_edge).second;
					//	
					//}
					//selected_edge = this->g.getEdge(selected_edge).second;
					//break;
				}
			}*/
		//}
		nextKey = this->g.getEdge(selected_edge).second;
		path.edgeKeys.push_back(selected_edge);
	}
	path.nodeKeys.push_back(nextKey);
	//for (std::vector<int>::iterator it1 = MCPath.nodeKeys.begin(); it1 != MCPath.nodeKeys.end(); it1++)
	//	std::cout << *it1 << std::endl;
	return true;
}

Path_t ConvexRelaxationSolver::getGreedyPath()
{	
	int nextKey = targetKey;
	Path_t greedyPath;
	while (nextKey != startKey)
	{
		greedyPath.nodeKeys.insert(greedyPath.nodeKeys.begin(),nextKey);
		std::vector<int> inEdges = this->g.findInEdges(nextKey);
		double max_y = 0.;
		int best_edge = -1;
		for (std::vector<int>::iterator it = inEdges.begin(); it != inEdges.end(); it++)
		{
			double current_y = this->relaxedSolution.y(*it);

			if (current_y > max_y)
			{
				int current_key = this->g.getEdge(*it).first;
				bool not_found = true;
				for (std::vector<int>::iterator itNodeKey = greedyPath.nodeKeys.begin(); itNodeKey != greedyPath.nodeKeys.end(); itNodeKey++)
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
		nextKey = this->g.getEdge(best_edge).first;
		greedyPath.edgeKeys.insert(greedyPath.edgeKeys.begin(),best_edge);
	}
	greedyPath.nodeKeys.insert(greedyPath.nodeKeys.begin(),nextKey);

	/*int nextKey = startKey;
	Path_t greedyPath;
	while (nextKey != targetKey)
	{
		greedyPath.nodeKeys.push_back(nextKey);
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
				for (std::vector<int>::iterator itNodeKey = greedyPath.nodeKeys.begin(); itNodeKey != greedyPath.nodeKeys.end(); itNodeKey++)
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
		greedyPath.edgeKeys.push_back(best_edge);
	}
	greedyPath.nodeKeys.push_back(nextKey);
	*/
	return greedyPath;
}

void ConvexRelaxationSolver::solve()
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
		this->relaxedSolution.x = Eigen::MatrixXd::Zero(this->nN, this->n);
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
}