#include<Eigen/Dense>
#include <chrono>
#include "Range.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "GCS.h"
#include "NavGraph.h"
#include "CObsConic.h"
#include "IRISConic.h"
#include "AStarIRISConic.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "ConvexRelaxationMinDistanceSolver.h"
#include "MinDistanceSolver.h"
#include <fstream>

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
	Eigen::Matrix<double, 2, 3> v3({ {0.,1.,2.},{-8.,-6.,-8.} });
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
	Eigen::Matrix<double, 2, 3> v3({ {0.,1.,2.},{-8.,-6.,-8.} });
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

int AStart_convex_relaxation_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar convex relaxation test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.addConvexSets(qstart);
	std::cout << "NavGraph of convex sets" << std::endl;
	irisConic.gcs.print();
	std::vector<int> edgeNodeKeys = irisConic.navGraph.getNodeKeys();
	for (int i = 0; i < irisConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = irisConic.navGraph.getNode(edgeNodeKeys[i]);
		if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PointNode*>(node) != NULL) {
			PointNode* pointNode = (PointNode*)node->getNodeData();
			pointNode->point.print();
		}
	}
	std::cout << "Edge graph of convex sets" << std::endl;
	irisConic.navGraph.print();
	ConvexRelaxationMinDistanceSolver solver(irisConic.navGraph,irisConic.qStartNodeNavGraphKey,irisConic.qTargetNodeNavGraphKey);
	solver.setTask();
	solver.solve();
	std::cout << "Relaxed Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.relaxedSolution.x << std::endl;
	std::cout << "y variables " << std::endl; 
	std::cout << solver.relaxedSolution.y << std::endl;
	std::cout << "z variables " << std::endl;
	std::cout << solver.relaxedSolution.z << std::endl;
	std::cout << "p variables " << std::endl;
	std::cout << solver.relaxedSolution.p << std::endl;
	std::cout << "l variables " << std::endl;
	std::cout << solver.relaxedSolution.l << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << solver.relaxedSolution.cost << std::endl;
	std::cout << "Feasible Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.feasibleSolution.x << std::endl;
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.nodeKeys.begin(); it != solver.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.edgeKeys.begin(); it != solver.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl; 
	std::cout << "cost " << std::endl;
	std::cout << solver.feasibleSolution.cost << std::endl;
	
	return 0;
}

int AStart_expand_convex_relaxation_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar expand convex relaxation test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.addConvexSets(qstart);
	std::cout << "Graph of convex sets" << std::endl;
	std::vector<int> nodeKeys = irisConic.gcs.getNodeKeys();
	for (int i = 0; i < irisConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = irisConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	irisConic.gcs.print();
	std::cout << "NavGraph of convex sets" << std::endl;
	std::vector<int> edgeNodeKeys = irisConic.navGraph.getNodeKeys();
	for (int i = 0; i < irisConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = irisConic.navGraph.getNode(edgeNodeKeys[i]);
		if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PointNode*>(node) != NULL) {
			PointNode* pointNode = (PointNode*)node->getNodeData();
			pointNode->point.print();
		}
	}
	irisConic.navGraph.print();
	irisConic.expandTerminalNode(2);
	std::cout << "Graph of convex sets" << std::endl;
	nodeKeys = irisConic.gcs.getNodeKeys();
	for (int i = 0; i < irisConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = irisConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	irisConic.gcs.print();
	std::cout << "NavGraph of convex sets" << std::endl;
	edgeNodeKeys = irisConic.navGraph.getNodeKeys();
	for (int i = 0; i < irisConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = irisConic.navGraph.getNode(edgeNodeKeys[i]);
		if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PointNode*>(node) != NULL) {
			PointNode* pointNode = (PointNode*)node->getNodeData();
			pointNode->point.print();
		}
	}
	irisConic.navGraph.print();
	ConvexRelaxationMinDistanceSolver solver(irisConic.navGraph, irisConic.qStartNodeNavGraphKey, irisConic.qTargetNodeNavGraphKey);
	solver.setTask();
	solver.solve();
	std::cout << "Relaxed Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.relaxedSolution.x << std::endl;
	std::cout << "y variables " << std::endl;
	std::cout << solver.relaxedSolution.y << std::endl;
	std::cout << "z variables " << std::endl;
	std::cout << solver.relaxedSolution.z << std::endl;
	std::cout << "p variables " << std::endl;
	std::cout << solver.relaxedSolution.p << std::endl;
	std::cout << "l variables " << std::endl;
	std::cout << solver.relaxedSolution.l << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << solver.relaxedSolution.cost << std::endl;
	std::cout << "Feasible Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.feasibleSolution.x << std::endl;
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.nodeKeys.begin(); it != solver.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.edgeKeys.begin(); it != solver.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << solver.feasibleSolution.cost << std::endl;

	return 0;
}

int AStart_expand_convex_MIP_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar expand convex MIP test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.addConvexSets(qstart);
	std::cout << "Graph of convex sets" << std::endl;
	std::vector<int> nodeKeys = irisConic.gcs.getNodeKeys();
	for (int i = 0; i < irisConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = irisConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	irisConic.gcs.print();
	std::cout << "NavGraph of convex sets" << std::endl;
	std::vector<int> edgeNodeKeys = irisConic.navGraph.getNodeKeys();
	for (int i = 0; i < irisConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = irisConic.navGraph.getNode(edgeNodeKeys[i]);
		if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PointNode*>(node) != NULL) {
			PointNode* pointNode = (PointNode*)node->getNodeData();
			pointNode->point.print();
		}
	}
	irisConic.navGraph.print();
	irisConic.expandTerminalNode(2);
	std::cout << "Graph of convex sets" << std::endl;
	nodeKeys = irisConic.gcs.getNodeKeys();
	for (int i = 0; i < irisConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = irisConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	irisConic.gcs.print();
	std::cout << "NavGraph of convex sets" << std::endl;
	edgeNodeKeys = irisConic.navGraph.getNodeKeys();
	for (int i = 0; i < irisConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = irisConic.navGraph.getNode(edgeNodeKeys[i]);
		if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PointNode*>(node) != NULL) {
			PointNode* pointNode = (PointNode*)node->getNodeData();
			pointNode->point.print();
		}
	}
	irisConic.navGraph.print();
	MinDistanceSolver solver(irisConic.navGraph, irisConic.qStartNodeNavGraphKey, irisConic.qTargetNodeNavGraphKey);
	solver.setTask();
	solver.solve();
	std::cout << "MIP Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.MIPSolution.x << std::endl;
	std::cout << "y variables " << std::endl;
	std::cout << solver.MIPSolution.y << std::endl;
	std::cout << "z variables " << std::endl;
	std::cout << solver.MIPSolution.z << std::endl;
	std::cout << "p variables " << std::endl;
	std::cout << solver.MIPSolution.p << std::endl;
	std::cout << "l variables " << std::endl;
	std::cout << solver.MIPSolution.l << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << solver.MIPSolution.cost << std::endl;
	std::cout << "Feasible Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.feasibleSolution.x << std::endl;
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.nodeKeys.begin(); it != solver.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.edgeKeys.begin(); it != solver.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << solver.feasibleSolution.cost << std::endl;

	return 0;
}


int AStar_IRIS_relaxed_solver_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar Relaxed Solver test" << std::endl;
	Range* range = Range::getInstance();
	double lb = -10.;
	double ub = 10.;
	range->setRange(2, lb, ub);
	CObsConic cObs = getCObsPolyhedronV();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISParams.ExpandableIRISParams.IRISParams.seperatingHyperplaneAligned = true;
	AStarIRISParams.ExpandableIRISParams.maxItersOptimalPath = 200;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::ofstream fout("AStar_IRIS_relaxed_solver_test_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	cObs.print(fout);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::duration<double, std::milli> ms_double;
	double phase1Time = 0.0;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.do_RelaxedSolver(fout);
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Resulting navGraph2 GCS" << std::endl;
	for (std::map<std::pair<int, int>, int>::iterator it = irisConic.navGraph2gcs.begin(); it != irisConic.navGraph2gcs.end(); it++)
	{
		std::cout << it->first.first << " " << it->first.second << " " << it->second << std::endl;
	}
	ms_double = t2 - t1;
	phase1Time = ms_double.count();
	std::cout << "Performance" << std::endl;
	std::cout << phase1Time << "ms" << std::endl;
	fout << "writerObj = VideoWriter('AStar_IRIS_relaxed_solver_test_results.avi');" << std::endl;
	fout << "writerObj.FrameRate = 10;" << std::endl;
	fout << "open(writerObj);" << std::endl;
	fout << "for i=1:length(Frames)" << std::endl;
	fout << "  frame = Frames(i) ;" << std::endl;
	fout << "  writeVideo(writerObj, frame);" << std::endl;
	fout << "end" << std::endl;
	fout << "close(writerObj);" << std::endl;
	return 0;
}

int AStar_IRIS_relaxed_solver_test1(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar Relaxed Solver test1" << std::endl;
	Range* range = Range::getInstance();
	double lb = -20.;
	double ub = 20.;
	range->setRange(2, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_15_obstacles();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISParams.ExpandableIRISParams.IRISParams.seperatingHyperplaneAligned = true;
	AStarIRISParams.ExpandableIRISParams.maxItersOptimalPath = 200;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::ofstream fout("AStar_IRIS_relaxed_solver_test1_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	cObs.print(fout);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::duration<double, std::milli> ms_double;
	double phase1Time = 0.0;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.do_RelaxedSolver(fout);
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	phase1Time = ms_double.count();
	std::cout << "Performance" << std::endl;
	std::cout << phase1Time << "ms" << std::endl;
	fout << "writerObj = VideoWriter('AStar_IRIS_relaxed_solver_test1_results.avi');" << std::endl;
	fout << "writerObj.FrameRate = 2;" << std::endl;
	fout << "open(writerObj);" << std::endl;
	fout << "for i=1:length(Frames)" << std::endl;
	fout << "  frame = Frames(i) ;" << std::endl;
	fout << "  writeVideo(writerObj, frame);" << std::endl;
	fout << "end" << std::endl;
	fout << "close(writerObj);" << std::endl;
	return 0;
}

int AStar_IRIS_relaxed_solver_test2(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar Relaxed Solver test2" << std::endl;
	Range* range = Range::getInstance();
	double lb = -35.;
	double ub = 35.;
	range->setRange(2, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_30_obstacles();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISParams.ExpandableIRISParams.IRISParams.seperatingHyperplaneAligned = true;
	AStarIRISParams.ExpandableIRISParams.maxItersOptimalPath = 300;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::ofstream fout("AStar_IRIS_relaxed_solver_test2_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	cObs.print(fout);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::duration<double, std::milli> ms_double;
	double phase1Time = 0.0;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.do_RelaxedSolver(fout);
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	phase1Time = ms_double.count();
	std::cout << "Performance" << std::endl;
	std::cout << phase1Time << "ms" << std::endl;
	fout << "writerObj = VideoWriter('AStar_IRIS_relaxed_solver_test2_results.avi');" << std::endl;
	fout << "writerObj.FrameRate = 2;" << std::endl;
	fout << "open(writerObj);" << std::endl;
	fout << "for i=1:length(Frames)" << std::endl;
	fout << "  frame = Frames(i) ;" << std::endl;
	fout << "  writeVideo(writerObj, frame);" << std::endl;
	fout << "end" << std::endl;
	fout << "close(writerObj);" << std::endl;
	return 0;
}

int AStar_IRIS_MIP_solver_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar MIP Solver test" << std::endl;
	Range* range = Range::getInstance();
	double lb = -10.;
	double ub = 10.;
	range->setRange(2, lb, ub);
	CObsConic cObs = getCObsPolyhedronV();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	//AStarIRISParams.ExpandableIRISParams.IRISParams.seperatingHyperplaneAligned = true;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::ofstream fout("AStar_IRIS_MIP_test_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	cObs.print(fout);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::duration<double, std::milli> ms_double;
	double phase1Time = 0.0;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.do_MIPSolver(fout);
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	phase1Time = ms_double.count();
	std::cout << "Performance" << std::endl;
	std::cout << phase1Time << "ms" << std::endl;
	/*std::cout << "Graph of convex sets" << std::endl;
	std::vector<int> nodeKeys = irisConic.gcs.getNodeKeys();
	for (int i = 0; i < irisConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = irisConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	irisConic.gcs.print();
	std::cout << "NavGraph of convex sets" << std::endl;
	std::vector<int> edgeNodeKeys = irisConic.navGraph.getNodeKeys();
	for (int i = 0; i < irisConic.navGraph.numNodes; i++)
	{
		Node* node = irisConic.navGraph.getNode(edgeNodeKeys[i]);
		if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL) {
			std::cout << "Terminal Node Convex set " << edgeNodeKeys[i] << std::endl;
			PolyhedronTerminalNode* polyNode = (PolyhedronTerminalNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
			std::cout << "Polyhedron Node Convex set " << edgeNodeKeys[i] << std::endl;
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PointNode*>(node) != NULL) {
			std::cout << "Point Node Convex set " << edgeNodeKeys[i] << std::endl;
			PointNode* pointNode = (PointNode*)node->getNodeData();
			pointNode->point.print();
		}
	}
	irisConic.navGraph.print();
	MinDistanceSolver solver(irisConic.navGraph, irisConic.qStartNodeNavGraphKey, irisConic.qTargetNodeNavGraphKey);
	solver.setTask();
	solver.solve();
	std::cout << "MIP Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.MIPSolution.x << std::endl;
	std::cout << "y variables " << std::endl;
	std::cout << solver.MIPSolution.y << std::endl;
	std::cout << "z variables " << std::endl;
	std::cout << solver.MIPSolution.z << std::endl;
	std::cout << "p variables " << std::endl;
	std::cout << solver.MIPSolution.p << std::endl;
	std::cout << "l variables " << std::endl;
	std::cout << solver.MIPSolution.l << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << solver.MIPSolution.cost << std::endl;
	std::cout << "Feasible Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.feasibleSolution.x << std::endl;
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.nodeKeys.begin(); it != solver.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.edgeKeys.begin(); it != solver.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << solver.feasibleSolution.cost << std::endl;*/

	return 0;
}

int main(int argc, char** argv)
{
	//terminalNode_randomseed_test();
	//AStart_qstart_qtarget_test(Eigen::Vector<double, 2>({ 0.,0. }), Eigen::Vector<double, 2>({9.,8}));
	//addConvexSets_test(Eigen::Vector<double, 2>({ 0.,0. }), Eigen::Vector<double, 2>({ 9.,8 }));
	//AStart_expand_terminal_node_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
	//AStart_convex_relaxation_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
	//AStart_expand_convex_relaxation_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
	//AStart_expand_convex_MIP_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
	//AStar_IRIS_MIP_solver_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
	AStar_IRIS_relaxed_solver_test(Eigen::Vector<double, 2>({ -5.,-8.}), Eigen::Vector<double, 2>({ 9.,8.}));
	//AStar_IRIS_relaxed_solver_test1(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8. }));
	//AStar_IRIS_relaxed_solver_test2(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8. }));
}