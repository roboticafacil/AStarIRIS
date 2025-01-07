#ifndef EXPANDABLE_IRIS_CONIC
#define EXPANDABLE_IRIS_CONIC
#include "Graph.h"
#include<Eigen/Dense>
#include "CObsConic.h"
#include "GCS.h"
#include "IRISConic.h"
#include "PointNode.h"
#include "ConvexRelaxationSolver.h"

typedef struct
{
	IRISParams_t IRISParams;
	int maxTrialsTerminal;
	int maxItersOptimalPath;
}ExpandableIRISParams_t;

class ExpandableIRISConic : public IRISConic
{
public:
	EGCS navGraph;
	int qStartNodeNavGraphKey;
	int qTargetNodeNavGraphKey;
public:
	//ExpandableIRISConic(const Graph& g, const std::vector<int>& nodeKeys, const int& startKey, const int& targetKey);
	ExpandableIRISConic(CObsConic& cObs, const ExpandableIRISParams_t& params);
	~ExpandableIRISConic();
	void addConvexSets(const Eigen::VectorXd& q);
	void buildNavGraph(const Eigen::VectorXd& qstart, const Eigen::VectorXd& qtarget);
	void terminalCosts(Eigen::VectorXd& hCosts, Eigen::MatrixXd& pClosest);
	void expandTerminalNode(const int& terminalNodeKey);
	Eigen::VectorXd generateRandomSeedInTerminalNode();
	static ExpandableIRISParams_t getDefaultExpandableIRISParams();
	EGCS getGraphWithoutTerminalConnections();
	void removeNode(const int key);
protected:
	std::map<std::pair<int, int>, int> egcs2gcs;
	std::vector<std::pair<int, int>> terminalEGCSNodeKeys;
	PointNode* qStartNode;
	PointNode* qTargetNode;
	ExpandableIRISParams_t params;
};
#endif