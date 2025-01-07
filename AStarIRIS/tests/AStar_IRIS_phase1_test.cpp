#include<Eigen/Dense>
#include <chrono>
#include "Range.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "GCS.h"
#include "EGCS.h"
#include "CObsConic.h"
#include "IRISConic.h"
#include "AStarIRISConic.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "ConvexRelaxationMinDistanceSolver.h"
#include "MinDistanceSolver.h"

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


int AStar_IRIS_phase1_relaxed_solver_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar phase1 test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::duration<double, std::milli> ms_double;
	double phase1Time = 0.0;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.doPhase1_RelaxedSolver();
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	phase1Time = ms_double.count();
	std::cout << "Phase 1 performance" << std::endl;
	std::cout << phase1Time << "ms" << std::endl;
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
	EGCS simplifiedGraph=irisConic.getGraphWithoutTerminalConnections();
	std::cout << "NavGraph of convex sets" << std::endl;
	std::vector<int> edgeNodeKeys = simplifiedGraph.getNodeKeys();
	for (int i = 0; i < simplifiedGraph.numNodes; i++)
	{
		
		Node* node = simplifiedGraph.getNode(edgeNodeKeys[i]);
		if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL) {
			std::cout << "Terminal Node Convex set " << edgeNodeKeys[i] << std::endl;
			PolyhedronTerminalNode* polyNode = (PolyhedronTerminalNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
			std::cout << "Intersection Convex set " << edgeNodeKeys[i] << std::endl;
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PointNode*>(node) != NULL) {
			std::cout << "Point Convex set " << edgeNodeKeys[i] << std::endl;
			PointNode* pointNode = (PointNode*)node->getNodeData();
			pointNode->point.print();
		}
	}
	simplifiedGraph.print();
	ConvexRelaxationMinDistanceSolver solver(simplifiedGraph, irisConic.qStartNodeNavGraphKey, irisConic.qTargetNodeNavGraphKey);
	solver.setTask();
	solver.solve();
	solver.computeFeasibleSolution(100);
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.nodeKeys.begin(); it != solver.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = solver.optimalPath.edgeKeys.begin(); it != solver.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;

	//Now compute a feasible solution using the MIP problem
	/*MinDistanceSolver solver(simplifiedGraph, irisConic.qStartNodeNavGraphKey, irisConic.qTargetNodeNavGraphKey);
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
	ConvexRelaxationMinDistanceSolver solver(simplifiedGraph, irisConic.qStartNodeNavGraphKey, irisConic.qTargetNodeNavGraphKey);
	solver.setTask();
	solver.solve();*/
	std::cout << "Feasible Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << solver.feasibleSolution.x << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << solver.feasibleSolution.cost << std::endl;
	return 0;
}

int AStar_IRIS_phase1_MIP_solver_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "AStar phase1 test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::duration<double, std::milli> ms_double;
	double phase1Time = 0.0;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.doPhase1_MIPSolver();
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	phase1Time = ms_double.count();
	std::cout << "Phase 1 performance" << std::endl;
	std::cout << phase1Time << "ms" << std::endl;
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
	std::cout << solver.feasibleSolution.cost << std::endl;

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
	AStar_IRIS_phase1_relaxed_solver_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
	//AStar_IRIS_phase1_MIP_solver_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
}