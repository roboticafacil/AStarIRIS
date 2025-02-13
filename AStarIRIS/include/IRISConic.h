#ifndef IRIS_CONIC
#define IRIS_CONIC
#include "Range.h"
#include "Node.h"
#include "CObsConic.h"
#include "GCS.h"
#include "NavGraph.h"

typedef struct
{
	int status;
}IRISRes_t;

/*typedef struct
{
	int nodeKey;
	int idx;
	Eigen::VectorXd ai;
	double bi;
}IRISNeighbour_t;*/

typedef struct
{
	int n;
	double tolConvexSetConvergence;
	double shrinkFactor;
	int maxTrials;
	bool seperatingHyperplaneAligned;
	double tol;
}IRISParams_t;

class IRISConic
{
public:
	GCS gcs;
	//IRISConic(CObsConic* cObs, GCS* gcs, EGCS* egcs, const IRISParams_t& params);
	IRISConic(CObsConic& cObs, const IRISParams_t& params);
	virtual std::vector<int> addConvexSets(const Eigen::VectorXd& q);
	virtual int addConvexSet(const Eigen::VectorXd& q);
	static IRISParams_t getDefaultIRISParams();
	//void IRISConic::separatingHyperplanes(Ellipsoid& ellipsoid, Polyhedron& polyhedron, std::vector<IRISNeighbour_t>& neighbours);
	void IRISConic::separatingHyperplanes(Ellipsoid& ellipsoid, Polyhedron& polyhedron);
	void generateGCS();
	void generateGCS(std::ostream& out);
	Eigen::VectorXd generateRandomSeed(int& trials);
	//void addStartNode(const Eigen::VectorXd& q);
public:
	IRISParams_t params;
protected:
	//void computeConvexSet(Ellipsoid& ellipsoid, Polyhedron& convexSet, std::vector<IRISNeighbour_t>& neighbours);
	virtual void computeConvexSet(Ellipsoid& ellipsoid, Polyhedron& convexSet, std::vector<int>& neighbourKeys);
protected:
	CObsConic* _cObs;
	
};
#endif