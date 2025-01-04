#ifndef A_STAR_IRIS_CONIC
#define A_STAR_IRIS_CONIC
#include "Graph.h"
#include<Eigen/Dense>
#include "CObsConic.h"
#include "GCS.h"
#include "IRISConic.h"
#include "PointNode.h"

typedef struct
{
	IRISParams_t IRISParams;
	int maxTrialsTerminal;
}AStarIRISParams_t;

class AStarIRISConic: public IRISConic
{
public:
	//EGCS egcs;
	EGCS navGraph;
	int qStartNodeNavGraphKey;
	int qTargetNodeNavGraphKey;
public:
	//AStarIRISConic(const Graph& g, const std::vector<int>& nodeKeys, const int& startKey, const int& targetKey);
	AStarIRISConic(CObsConic& cObs, const AStarIRISParams_t& params);
	void addConvexSets(const Eigen::VectorXd& q);
	void buildNavGraph(const Eigen::VectorXd& qstart, const Eigen::VectorXd& qtarget);
	void terminalCosts(Eigen::VectorXd& hCosts, Eigen::MatrixXd& pClosest);
	void expandTerminalNode(const int &terminalNodeKey);
	Eigen::VectorXd generateRandomSeedInTerminalNode();
	static AStarIRISParams_t getDefaultAStarIRISParams();
	void doPhase1();
private:
	std::map<std::pair<int, int>, int> egcs2gcs;
	std::vector<std::pair<int, int>> terminalEGCSNodeKeys;
	PointNode* qStartNode;
	PointNode* qTargetNode;
	//Graph graph;
	//std::vector<int> nodeKeys;
	AStarIRISParams_t params;
};
#endif