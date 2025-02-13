#include<Eigen/Dense>
#include <chrono>
#include "Range.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "GCS.h"
#include "CObsConic.h"
#include "IRISConic.h"
#include "PolyhedronNode.h"
#include <fstream>

CObsConic getCObsPolyhedronV()
{
	CObsConic cObs = CObsConic();
	Eigen::Matrix<double, 2, 4> v0({ {-3.,-4.,-4.,-5.},{5.,6.,4.,5.} });
	PolyhedronV* o0 = new PolyhedronV(v0);
	cObs.addObject(o0);
	Eigen::Matrix<double, 2, 5> v1({ {5.,6.,8.,9.,8.},{4.,8.,7.,5.,2.} });
	PolyhedronV* o1 = new PolyhedronV(v1);
	cObs.addObject(o1);
	Eigen::Matrix<double, 2, 4> v2({ {-9.,-7.,-1.,-2.},{-7.,-2.,-3.,-4.} });
	PolyhedronV* o2 = new PolyhedronV(v2);
	cObs.addObject(o2);
	Eigen::Matrix<double, 2, 3> v3({ {0.,1.,2.},{-9.,-7.,-9.} });
	PolyhedronV* o3 = new PolyhedronV(v3);
	cObs.addObject(o3);
	Eigen::Matrix<double, 2, 5> v4({ {5.,5.,6.,7.,7.},{-7.,-5.,-4.,-5.,-7.} });
	PolyhedronV* o4 = new PolyhedronV(v4);
	cObs.addObject(o4);
	return cObs;
}

CObsConic getCObsPolyhedronV_15_obstacles()
{
	CObsConic cObs = CObsConic();
	Eigen::Matrix<double, 2, 4> v0({ {-3.,-4.,-4.,-5.},{5.,6.,4.,5.} });
	PolyhedronV* o0 = new PolyhedronV(v0);
	cObs.addObject(o0);
	Eigen::Matrix<double, 2, 5> v1({ {5.,6.,8.,9.,8.},{4.,8.,7.,5.,2.} });
	PolyhedronV* o1 = new PolyhedronV(v1);
	cObs.addObject(o1);
	Eigen::Matrix<double, 2, 4> v2({ {-9.,-7.,-1.,-2.},{-7.,-2.,-3.,-4.} });
	PolyhedronV* o2 = new PolyhedronV(v2);
	cObs.addObject(o2);
	Eigen::Matrix<double, 2, 3> v3({ {0.,1.,2.},{-9.,-7.,-9.} });
	PolyhedronV* o3 = new PolyhedronV(v3);
	cObs.addObject(o3);
	Eigen::Matrix<double, 2, 5> v4({ {5.,5.,6.,7.,7.},{-7.,-5.,-4.,-5.,-7.} });
	PolyhedronV* o4 = new PolyhedronV(v4);
	cObs.addObject(o4);
	Eigen::Matrix<double, 2, 5> v5({ {-20., -18., -10., -8., -10},{-5., -3., -18., -16., -10.} });
	PolyhedronV* o5 = new PolyhedronV(v5);
	cObs.addObject(o5);
	Eigen::Matrix<double, 2, 4> v6({ {0., 15., 15., 0.},{-16., -16., -13., -13.} });
	PolyhedronV* o6 = new PolyhedronV(v6);
	cObs.addObject(o6);
	Eigen::Matrix<double, 2, 4> v7({ {-6., -5., -5., -7.},{-11., -11., -10., -10.} });
	PolyhedronV* o7 = new PolyhedronV(v7);
	cObs.addObject(o7);
	Eigen::Matrix<double, 2, 4> v8({ {15., 16., 16., 15.},{-10., -10., -5., -5.} });
	PolyhedronV* o8 = new PolyhedronV(v8);
	cObs.addObject(o8);
	Eigen::Matrix<double, 2, 4> v9({ {12., 13., 13., 12.},{3., 5., 15., 15.} });
	PolyhedronV* o9 = new PolyhedronV(v9);
	cObs.addObject(o9);
	Eigen::Matrix<double, 2, 4> v10({ {-12., -13., -13., -12.},{3., 5., 15., 15.} });
	PolyhedronV* o10 = new PolyhedronV(v10);
	cObs.addObject(o10);
	Eigen::Matrix<double, 2, 4> v11({ {-5., 5., 5., -5.},{9., 9., 10., 10.} });
	PolyhedronV* o11 = new PolyhedronV(v11);
	cObs.addObject(o11);
	Eigen::Matrix<double, 2, 4> v12({ {9., 10., 10., 9.},{-5, -5., 2., 2.} });
	PolyhedronV* o12 = new PolyhedronV(v12);
	cObs.addObject(o12);
	Eigen::Matrix<double, 2, 4> v13({ {-7., -8., -8., -9.},{6., 7., 5., 6.} });
	PolyhedronV* o13 = new PolyhedronV(v13);
	cObs.addObject(o13);
	Eigen::Matrix<double, 2, 4> v14({ {-12., -13., -13., -14.},{-1., 0., -2., -1.} });
	PolyhedronV* o14 = new PolyhedronV(v14);
	cObs.addObject(o14);
	//Eigen::Matrix<double, 2, 4> v15({ {0., 8., 8., 0.},{12., 12., 13., 13.} });
	//PolyhedronV* o15 = new PolyhedronV(v15);
	//cObs.addObject(o15);
	return cObs;
}


CObsConic getCObsPolyhedronV_30_obstacles()
{
	CObsConic cObs = CObsConic();
	Eigen::Matrix<double, 2, 4> v0({ {-3.,-4.,-4.,-5.},{5.,6.,4.,5.} });
	PolyhedronV* o0 = new PolyhedronV(v0);
	cObs.addObject(o0);
	Eigen::Matrix<double, 2, 5> v1({ {5.,6.,8.,9.,8.},{4.,8.,7.,5.,2.} });
	PolyhedronV* o1 = new PolyhedronV(v1);
	cObs.addObject(o1);
	Eigen::Matrix<double, 2, 4> v2({ {-9.,-7.,-1.,-2.},{-7.,-2.,-3.,-4.} });
	PolyhedronV* o2 = new PolyhedronV(v2);
	cObs.addObject(o2);
	Eigen::Matrix<double, 2, 3> v3({ {0.,1.,2.},{-9.,-7.,-9.} });
	PolyhedronV* o3 = new PolyhedronV(v3);
	cObs.addObject(o3);
	Eigen::Matrix<double, 2, 5> v4({ {5.,5.,6.,7.,7.},{-7.,-5.,-4.,-5.,-7.} });
	PolyhedronV* o4 = new PolyhedronV(v4);
	cObs.addObject(o4);
	Eigen::Matrix<double, 2, 5> v5({ {-20., -18., -10., -8., -10},{-5., -3., -18., -16., -10.} });
	PolyhedronV* o5 = new PolyhedronV(v5);
	cObs.addObject(o5);
	Eigen::Matrix<double, 2, 4> v6({ {0., 15., 15., 0.},{-16., -16., -13., -13.} });
	PolyhedronV* o6 = new PolyhedronV(v6);
	cObs.addObject(o6);
	Eigen::Matrix<double, 2, 4> v7({ {-6., -5., -5., -7.},{-11., -11., -10., -10.} });
	PolyhedronV* o7 = new PolyhedronV(v7);
	cObs.addObject(o7);
	Eigen::Matrix<double, 2, 4> v8({ {15., 16., 16., 15.},{-10., -10., -5., -5.} });
	PolyhedronV* o8 = new PolyhedronV(v8);
	cObs.addObject(o8);
	Eigen::Matrix<double, 2, 4> v9({ {12., 13., 13., 12.},{3., 5., 15., 15.} });
	PolyhedronV* o9 = new PolyhedronV(v9);
	cObs.addObject(o9);
	Eigen::Matrix<double, 2, 4> v10({ {-12., -13., -13., -12.},{3., 5., 15., 15.} });
	PolyhedronV* o10 = new PolyhedronV(v10);
	cObs.addObject(o10);
	Eigen::Matrix<double, 2, 4> v11({ {-5., 5., 5., -5.},{9., 9., 10., 10.} });
	PolyhedronV* o11 = new PolyhedronV(v11);
	cObs.addObject(o11);
	Eigen::Matrix<double, 2, 4> v12({ {9., 10., 10., 9.},{-5, -5., 2., 2.} });
	PolyhedronV* o12 = new PolyhedronV(v12);
	cObs.addObject(o12);
	Eigen::Matrix<double, 2, 4> v13({ {-7., -8., -8., -9.},{6., 7., 5., 6.} });
	PolyhedronV* o13 = new PolyhedronV(v13);
	cObs.addObject(o13);
	Eigen::Matrix<double, 2, 4> v14({ {-12., -13., -13., -14.},{-1., 0., -2., -1.} });
	PolyhedronV* o14 = new PolyhedronV(v14);
	cObs.addObject(o14);
	Eigen::Matrix<double, 2, 4> v15({ {-8., 0., 0., -8.},{12., 12., 13., 13.} });
	PolyhedronV* o15 = new PolyhedronV(v15);
	cObs.addObject(o15);
	Eigen::Matrix<double, 2, 3> v16({ {-10., -5., 0.},{-20.,-25.,-20.} });
	PolyhedronV* o16 = new PolyhedronV(v16);
	cObs.addObject(o16);
	Eigen::Matrix<double, 2, 3> v17({ {2.,2.,12.},{20.,15.,20.} });
	PolyhedronV* o17 = new PolyhedronV(v17);
	cObs.addObject(o17);
	Eigen::Matrix<double, 2, 3> v18({ {15., 20., 25.},{-20.,-25.,-20.} });
	PolyhedronV* o18 = new PolyhedronV(v18);
	cObs.addObject(o18);
	Eigen::Matrix<double, 2, 3> v19({ {10., 5., 0.},{-30.,-25.,-30.} });
	PolyhedronV* o19 = new PolyhedronV(v19);
	cObs.addObject(o19);
	Eigen::Matrix<double, 2, 3> v20({ {-30., -30., -20.},{-20.,-30.,-30.} });
	PolyhedronV* o20 = new PolyhedronV(v20);
	cObs.addObject(o20);
	Eigen::Matrix<double, 2, 3> v21({ {-30., -30., -20.},{20.,30.,30.} });
	PolyhedronV* o21 = new PolyhedronV(v21);
	cObs.addObject(o21);
	Eigen::Matrix<double, 2, 3> v22({ {20., 30., 30.},{20.,20.,30.} });
	PolyhedronV* o22 = new PolyhedronV(v22);
	cObs.addObject(o22);
	Eigen::Matrix<double, 2, 4> v23({ {-23.,-24.,-24.,-25.},{15.,16.,14.,15.} });
	PolyhedronV* o23 = new PolyhedronV(v23);
	cObs.addObject(o23);
	Eigen::Matrix<double, 2, 4> v24({ {-23.,-24.,-24.,-25.},{-15.,-16.,-14.,-15.} });
	PolyhedronV* o24 = new PolyhedronV(v24);
	cObs.addObject(o24);
	Eigen::Matrix<double, 2, 4> v25({ {23.,24.,24.,25.},{0.,-1.,1.,0.} });
	PolyhedronV* o25 = new PolyhedronV(v25);
	cObs.addObject(o25);
	Eigen::Matrix<double, 2, 4> v26({ {17.,18.,18.,19.},{-15.,-16.,-14.,-15.} });
	PolyhedronV* o26 = new PolyhedronV(v26);
	cObs.addObject(o26);
	Eigen::Matrix<double, 2, 4> v27({ {15.,16.,16.,17.},{15.,16.,14.,15.} });
	PolyhedronV* o27 = new PolyhedronV(v27);
	cObs.addObject(o27);
	Eigen::Matrix<double, 2, 5> v28({ {-15.,-16.,-18.,-19.,-18.},{4.,8.,7.,5.,2.} });
	PolyhedronV* o28 = new PolyhedronV(v28);
	cObs.addObject(o28);
	Eigen::Matrix<double, 2, 5> v29({ {-5.,-6.,-8.,-9.,-8.},{24.,28.,27.,25.,22.} });
	PolyhedronV* o29 = new PolyhedronV(v29);
	cObs.addObject(o29);
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
	Eigen::Matrix<double, 2, 3> v3({ {0.,1.,2.},{-9.,-7.,-9.} });
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

int IRIS_test1(const double& shrinkFactor, const int& maxTrials, const double& tol)
{
	std::cout << "Generate GCS PolyhedronV test" << std::endl;
	Range* range = Range::getInstance();
	double lb = -20.;
	double ub = 20.;
	range->setRange(2, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_15_obstacles();
	IRISParams_t& IRISParams = IRISConic::getDefaultIRISParams();
	IRISParams.n = 2;
	IRISParams.shrinkFactor = shrinkFactor;
	IRISParams.maxTrials = maxTrials;
	IRISParams.tol = tol;
	//IRISParams.seperatingHyperplaneAligned = true;
	std::ofstream fout("IRIS_generation_test1_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	cObs.print(fout);
	IRISConic irisConic = IRISConic(cObs, IRISParams);
	std::chrono::duration<double, std::milli> ms_double;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.generateGCS(fout);
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
	fout << "writerObj = VideoWriter('IRIS_generation_test1_results.avi');" << std::endl;
	fout << "writerObj.FrameRate = 2;" << std::endl;
	fout << "open(writerObj);" << std::endl;
	fout << "for i=1:length(Frames)" << std::endl;
	fout << "  frame = Frames(i) ;" << std::endl;
	fout << "  writeVideo(writerObj, frame);" << std::endl;
	fout << "end" << std::endl;
	fout << "close(writerObj);" << std::endl;
	return 0;
}

int main(int argc, char** argv)
{
	//addConvexSets_test(Eigen::Vector<double, 2>({0.,0.}));
	//addConvexSets_test(Eigen::Vector<double, 2>({5.,-4.9}));
	//addConvexSets_test(Eigen::Vector<double, 2>({0.,0.}), Eigen::Vector<double, 2>({ 4.,-8. }), Eigen::Vector<double, 2>({4.5,-5.}));
	//addConvexSets_test(Eigen::Vector<double, 2>({ 5.,-4.9 }), Eigen::Vector<double, 2>({4.,-8.}), Eigen::Vector<double, 2>({ 4.5,-6.}));
	//addConvexSetsCircularRobot_test(Eigen::Vector<double, 2>({ 0.,0. })); //Este todavía no funciona bien!!
	IRIS_test1(1.,999,1e-3);
}