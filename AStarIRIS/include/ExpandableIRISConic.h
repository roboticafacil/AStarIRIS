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
	NavGraph navGraph;
	int qStartNodeNavGraphKey;
	int qTargetNodeNavGraphKey;
	std::map<std::pair<int, int>, int> navGraph2gcs;
public:
	//ExpandableIRISConic(const Graph& g, const std::vector<int>& nodeKeys, const int& startKey, const int& targetKey);
	ExpandableIRISConic(CObsConic& cObs, const ExpandableIRISParams_t& params);
	~ExpandableIRISConic();
	virtual int addConvexSet(const Eigen::VectorXd& q);
	virtual std::vector<int> addConvexSets(const Eigen::VectorXd& q);
	void buildNavGraph(const Eigen::VectorXd& qstart, const Eigen::VectorXd& qtarget);
	void terminalCosts(Eigen::VectorXd& hCosts, Eigen::MatrixXd& pClosest);
	void expandTerminalNode(const int& terminalNodeKey);
	Eigen::VectorXd generateRandomSeedInTerminalNode();
	static ExpandableIRISParams_t getDefaultExpandableIRISParams();
	NavGraph getGraphWithoutTerminalConnections();
	void removeNode(const int key);
protected:
	
	std::vector<std::pair<int, int>> terminalNodeKeys;
	PointNode* qStartNode;
	PointNode* qTargetNode;
	ExpandableIRISParams_t params;
private:
	void computeConvexSet(Ellipsoid& ellipsoid, Polyhedron& convexSet, std::vector<int>& neighbourKeys);
};
#endif