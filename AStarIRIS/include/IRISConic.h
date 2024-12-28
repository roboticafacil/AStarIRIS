#ifndef IRIS_CONIC
#define IRIS_CONIC
#include "Range.h"
#include "Node.h"
#include "CObsConic.h"
#include "GCS.h"

typedef struct
{
	int status;
}IRISRes_t;

typedef struct
{
	int nodeKey;
	Eigen::VectorXd ai;
	double bi;
}IRISNeighbour_t;

typedef struct
{
	int n;
	double tolConvexSetConvergence;
	double shrinkFactor;
	int maxCount;
}IRISParams_t;

class IRISConic
{
public:
	GCS* gcs;
	CObsConic* cObs;
	IRISConic(CObsConic* cObs, GCS* gcs, const IRISParams_t& params);
	void addConvexSets(const Eigen::VectorXd& q);
	static IRISParams_t getDefaultParams();
	void IRISConic::separatingHyperplanes(Ellipsoid& ellipsoid, Polyhedron* polyhedron, std::vector<IRISNeighbour_t>& neighbours);
	void generateGCS();
	Eigen::VectorXd generateRandomSeed(int& count);
public:
	IRISParams_t params;
private:
	
};
#endif