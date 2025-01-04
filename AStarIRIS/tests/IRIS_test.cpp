#include<Eigen/Dense>
#include <chrono>
#include "Range.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "GCS.h"
#include "CObsConic.h"
#include "IRISConic.h"
#include "PolyhedronNode.h"

CObsConic getCObsPolyhedronV()
{
	CObsConic cObs = CObsConic();
	Eigen::Matrix<double, 2, 4> v0({{-3.,-4.,-4.,-5.},{5.,6.,4.,5.}});
	PolyhedronV* o0 = new PolyhedronV(v0);
	cObs.addObject(o0);
	Eigen::Matrix<double, 2, 5> v1({{5.,6.,8.,9.,8.},{4.,8.,7.,5.,2.}});
	PolyhedronV* o1 = new PolyhedronV(v1);
	cObs.addObject(o1);
	Eigen::Matrix<double, 2, 4> v2({{-9.,-7.,-1.,-2.},{-7.,-2.,-3.,-4.}});
	PolyhedronV* o2 = new PolyhedronV(v2);
	cObs.addObject(o2);
	Eigen::Matrix<double, 2, 3> v3({{0.,1.,2.},{-8.,-6.,-8.}});
	PolyhedronV* o3 = new PolyhedronV(v3);
	cObs.addObject(o3);
	Eigen::Matrix<double, 2, 5> v4({{5.,5.,6.,7.,7.},{-7.,-5.,-4.,-5.,-7.}});
	PolyhedronV* o4 = new PolyhedronV(v4);
	cObs.addObject(o4);
	return cObs;
}

CObsConic getCObsCircularRobot(const double &radius)
{
	CObsConic cObs = CObsConic();
	Eigen::Matrix<double, 2, 4> v0({ {-3.,-4.,-4.,-5.},{5.,6.,4.,5.} });
	PolyhedronObstacleCircularRobot* o0 = new PolyhedronObstacleCircularRobot(v0,radius);
	cObs.addObject(o0);
	Eigen::Matrix<double, 2, 5> v1({ {5.,6.,8.,9.,8.},{4.,8.,7.,5.,2.} });
	PolyhedronObstacleCircularRobot* o1 = new PolyhedronObstacleCircularRobot(v1,radius);
	cObs.addObject(o1);
	Eigen::Matrix<double, 2, 4> v2({ {-9.,-7.,-1.,-2.},{-7.,-2.,-3.,-4.} });
	PolyhedronObstacleCircularRobot* o2 = new PolyhedronObstacleCircularRobot(v2,radius);
	cObs.addObject(o2);
	Eigen::Matrix<double, 2, 3> v3({ {0.,1.,2.},{-8.,-6.,-8.} });
	PolyhedronObstacleCircularRobot* o3 = new PolyhedronObstacleCircularRobot(v3,radius);
	cObs.addObject(o3);
	Eigen::Matrix<double, 2, 5> v4({ {5.,5.,6.,7.,7.},{-7.,-5.,-4.,-5.,-7.} });
	PolyhedronObstacleCircularRobot* o4 = new PolyhedronObstacleCircularRobot(v4,radius);
	cObs.addObject(o4);
	return cObs;
}

int addConvexSets_test(Eigen::Vector<double, 2>& q)
{
	std::cout << "Add ConvexSets PolyhedronV test" << std::endl;
	std::cout << "Testing seed " << std::endl;
	std::cout << q << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	IRISParams_t& IRISParams = IRISConic::getDefaultIRISParams();
	IRISParams.n = 2;
	IRISConic irisConic = IRISConic(cObs, IRISParams);
	irisConic.addConvexSets(q);
	std::vector<int> nodeKeys = irisConic.gcs.getNodeKeys();
	std::cout << "Results" << std::endl;
	for (int i = 0; i < irisConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = irisConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	std::cout << "Graph of convex sets" << std::endl;
	irisConic.gcs.print();
	return 0;
}


int addConvexSetsCircularRobot_test(Eigen::Vector<double, 2>& q)
{
	std::cout << "Add ConvexSets Circular Robot test" << std::endl;
	std::cout << "Testing seed " << std::endl;
	std::cout << q << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	double radius = 0.3;
	CObsConic cObs = getCObsCircularRobot(radius);
	IRISParams_t& IRISParams = IRISConic::getDefaultIRISParams();
	IRISParams.n = 2;
	IRISConic irisConic = IRISConic(cObs, IRISParams);
	irisConic.addConvexSets(q);
	std::vector<int> nodeKeys = irisConic.gcs.getNodeKeys();
	std::cout << "Results" << std::endl;
	for (int i = 0; i < irisConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = irisConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	std::cout << "Graph of convex sets" << std::endl;
	irisConic.gcs.print();
	return 0;
}

int IRIS_test(const double & shrinkFactor, const int &maxTrials, const double &tol)
{
	std::cout << "Generate GCS PolyhedronV test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	IRISParams_t& IRISParams = IRISConic::getDefaultIRISParams();
	IRISParams.n = 2;
	IRISParams.shrinkFactor = shrinkFactor;
	IRISParams.maxTrials = maxTrials;
	IRISParams.tol = tol;
	IRISConic irisConic = IRISConic(cObs, IRISParams);
	std::chrono::duration<double, std::milli> ms_double;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.generateGCS();
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	std::cout << "Total time " << ms_double.count() << " ms" << std::endl;
	std::vector<int> nodeKeys = irisConic.gcs.getNodeKeys();
	std::cout << "Results" << std::endl;
	for (int i = 0; i < irisConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = irisConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	std::cout << "Graph of convex sets" << std::endl;
	irisConic.gcs.print();
	return 0;
}

int main(int argc, char** argv)
{
	addConvexSets_test(Eigen::Vector<double, 2>({0.,0.}));
	addConvexSets_test(Eigen::Vector<double, 2>({5.,-4.9}));
	//addConvexSets_test(Eigen::Vector<double, 2>({0.,0.}), Eigen::Vector<double, 2>({ 4.,-8. }), Eigen::Vector<double, 2>({4.5,-5.}));
	//addConvexSets_test(Eigen::Vector<double, 2>({ 5.,-4.9 }), Eigen::Vector<double, 2>({4.,-8.}), Eigen::Vector<double, 2>({ 4.5,-6.}));
	//addConvexSetsCircularRobot_test(Eigen::Vector<double, 2>({ 0.,0. })); //Este todavía no funciona bien!!
	IRIS_test(1.,999,1e-3);
}