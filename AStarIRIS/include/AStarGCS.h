#ifndef A_STAR_GCS
#define A_STAR_GCS
#include "Graph.h"
#include<Eigen/Dense>

class AStarGCS
{
private:
	Eigen::VectorXd h;
	Eigen::MatrixXd pClosest;
	Graph graph;
	std::vector<int> nodeKeys;
	int startKey;
	int targetKey;
public:
	AStarGCS(const Graph& g, const std::vector<int>& nodeKeys, const int& startKey, const int& targetKey);
	Eigen::VectorXd getHeuristicCosts();
private:
	void computeHeuristics();
	int nN;
	int n;
};
#endif