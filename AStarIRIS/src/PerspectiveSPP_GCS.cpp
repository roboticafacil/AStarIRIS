#include "PerspectiveSPP_GCS.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PointNode.h"
#include "fusion.h"
#include <vector>
#include "EigenUtils.h"
#include <random>
#include <limits>

PerspectiveSPP_GCS::PerspectiveSPP_GCS(Graph* g, const int &N, const int& startKey, const int& targetKey) : SPP_GCS(g,N,startKey,targetKey), graphSimplified(false), perspectiveSolutionSolved(false)
{
}

bool PerspectiveSPP_GCS::getMCPath(Path_t & path)
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
		std::vector<int> outEdges = this->graphSimplified? this->simplifiedGraph.findOutEdges(nextKey) : this->g->findOutEdges(nextKey);
		std::vector<double> weights;
		for (std::vector<int>::iterator it = outEdges.begin(); it != outEdges.end(); it++)
		{
			weights.push_back(this->perspectiveSolution.y(*it));
		}
		std::discrete_distribution<int> dist(weights.begin(), weights.end());

		
		int selected_edge;
		//bool found = true;
		
		//while (found)
		//{
			std::mt19937 rng(rd());
			selected_edge = outEdges[dist(rng)];
			nextKey = this->graphSimplified ? this->simplifiedGraph.getEdge(selected_edge).second : this->g->getEdge(selected_edge).second;
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
		nextKey = this->graphSimplified ? this->simplifiedGraph.getEdge(selected_edge).second : this->g->getEdge(selected_edge).second;
		path.edgeKeys.push_back(selected_edge);
	}
	path.nodeKeys.push_back(nextKey);
	//for (std::vector<int>::iterator it1 = MCPath.nodeKeys.begin(); it1 != MCPath.nodeKeys.end(); it1++)
	//	std::cout << *it1 << std::endl;
	return true;
}

Path_t PerspectiveSPP_GCS::getGreedyPath()
{	
	int nextKey = targetKey;
	Path_t greedyPath;
	while (nextKey != startKey)
	{
		greedyPath.nodeKeys.insert(greedyPath.nodeKeys.begin(),nextKey);
		std::vector<int> inEdges = this->g->findInEdges(nextKey);
		double max_y = 0.;
		int best_edge = -1;
		for (std::vector<int>::iterator it = inEdges.begin(); it != inEdges.end(); it++)
		{
			double current_y = this->perspectiveSolution.y(*it);

			if (current_y > max_y)
			{
				int current_key = this->g->getEdge(*it).first;
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
		nextKey = this->g->getEdge(best_edge).first;
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

void PerspectiveSPP_GCS::simplifyGraph()
{
	if (this->perspectiveSolutionSolved)
	{
		//Create a copy of g without nodes and edges with very little probability of being selected based on the relaxed solution
		this->simplifiedGraph = Graph(g);
		std::vector<Edge> edges = this->g->getEdges();
		std::vector<double> simplifiedWeights;
		for (int i = 0; i < this->perspectiveSolution.y.rows(); i++)
		{
			if (this->perspectiveSolution.y(i) >= 1e-5)
			{
				simplifiedWeights.push_back(this->perspectiveSolution.y(i));
			}
			else
			{
				this->simplifiedGraph.removeEdge(edges[i].first, edges[i].second);
			}
		}
		std::vector<int> nodeKeys = this->simplifiedGraph.getNodeKeys();
		std::vector<int> nodeKeysToRemove;
		for (int i = 0; i < this->simplifiedGraph.numNodes; i++)
		{
			if (nodeKeys[i] != targetKey)
			{
				std::vector<int> outEdges = this->simplifiedGraph.findOutEdges(nodeKeys[i]);
				if (outEdges.size() == 0)
				{
					nodeKeysToRemove.push_back(nodeKeys[i]);
				}
			}
			else if (nodeKeys[i] != startKey)
			{
				std::vector<int> inEdges = this->simplifiedGraph.findInEdges(nodeKeys[i]);
				if (inEdges.size() == 0)
				{
					nodeKeysToRemove.push_back(nodeKeys[i]);
				}
			}
		}
		for (std::vector<int>::iterator it = nodeKeysToRemove.begin(); it != nodeKeysToRemove.end(); it++)
		{
			this->simplifiedGraph.removeNode(*it);
		}
		this->graphSimplified = true;
		this->simplifiedGraph.print();
	}
}