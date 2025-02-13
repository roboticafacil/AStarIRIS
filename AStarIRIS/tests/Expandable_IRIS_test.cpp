#include<Eigen/Dense>
#include <chrono>
#include "Range.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "GCS.h"
#include "NavGraph.h"
#include "CObsConic.h"
#include "IRISConic.h"
#include "ExpandableIRISConic.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "ConvexRelaxationMinDistanceSPP_GCS.h"
#include "MIPMinDistanceSPP_GCS.h"
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


int qstart_qtarget_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "Expandable qstart and qtarget test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	ExpandableIRISParams_t& ExpandableIRISParams = ExpandableIRISConic::getDefaultExpandableIRISParams();
	ExpandableIRISParams.IRISParams.n = 2;
	ExpandableIRISConic expandableIRISConic = ExpandableIRISConic(cObs, ExpandableIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	expandableIRISConic.buildNavGraph(qstart, qtarget);
	expandableIRISConic.addConvexSets(qstart);
	Eigen::VectorXd h;
	Eigen::MatrixXd pClosest;
	expandableIRISConic.terminalCosts(h, pClosest);
	std::vector<int> nodeKeys = expandableIRISConic.gcs.getNodeKeys();
	std::cout << "Results" << std::endl;
	for (int i = 0; i < expandableIRISConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = expandableIRISConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	std::cout << "NavGraph of convex sets" << std::endl;
	expandableIRISConic.gcs.print();
	std::vector<int> edgeNodeKeys = expandableIRISConic.navGraph.getNodeKeys();
	for (int i = 0; i < expandableIRISConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = expandableIRISConic.navGraph.getNode(edgeNodeKeys[i]);
		if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		else if (dynamic_cast<PointNode*>(node) != NULL) {
			PointNode* pointNode = (PointNode*)node->getNodeData();
			pointNode->point.print();
		}
	}
	std::cout << "Heuristic costs " << std::endl;
	std::cout << h << std::endl;
	std::cout << "Closest points in Terminal Nodes " << std::endl;
	std::cout << pClosest << std::endl;
	std::cout << "Edge graph of convex sets" << std::endl;
	expandableIRISConic.navGraph.print();
	return 0;
}

int terminalNode_randomseed_test()
{
	std::cout << "Expandable terminal node random seed test" << std::endl;
	int terminalConstrainIdx = 2;
	Eigen::Matrix<double, 5, 3> A({ {0.6996,-0.6487,0.2996},{0.5703,0.6724,0.4719},{0.,0.8723,-0.4889},{-0.5349,0.6181,-0.5760},{0.0480,0.9979,-0.0427} });
	Eigen::Vector<double, 5> b({ 3.5689,-1.3095,-3.4443,0.9384,-0.1748 });
	int n = A.cols();
	int m = A.rows();
	Eigen::VectorXd ai = A.row(terminalConstrainIdx);
	double bi = b(terminalConstrainIdx);
	int col;
	double aipivot = ai.cwiseAbs().maxCoeff(&col);
	std::vector<int> rows;
	for (int i = 0; i < m; i++)
	{
		if (i != terminalConstrainIdx)
			rows.push_back(i);
	}
	std::vector<int> cols;
	for (int i = 0; i < n; i++)
	{
		if (i != col)
			cols.push_back(i);
	}
	Eigen::MatrixXd Aj = A(rows, Eigen::placeholders::all);
	Eigen::VectorXd bj = b(rows);
	Eigen::VectorXd ai_col = ai(cols);
	std::cout << "A=" << A << std::endl;
	std::cout << "b=" << b << std::endl;
	std::cout << "ai=" << ai << std::endl;
	std::cout << "aipivot=" << aipivot << std::endl;
	std::cout << "ai_col=" << ai_col << std::endl;
	std::cout << "Aj=" << Aj << std::endl;
	std::cout << "bj=" << bj << std::endl;

	Range* range = Range::getInstance();
	std::pair<double, double> bounds = range->getBounds();
	Eigen::VectorXd seed;
	seed = Eigen::VectorXd::Random(n, 1);
	seed = (seed + Eigen::VectorXd::Constant(n, 1, 1.)) * (bounds.second - bounds.first) / 2.;
	seed = (seed + Eigen::VectorXd::Constant(n, 1, bounds.first));
	std::cout << "seed=" << seed << std::endl;
	seed(col) = (bi + (2 * 1e-3) - ai_col.dot(seed(cols))) / aipivot;
	std::cout << "projected seed=" << seed << std::endl;
	std::cout << "Contraint LHS " << (Aj * seed - bj) << std::endl;
	bool constrain = ((Aj * seed - bj).array() <= 0.).all();
	if (constrain)
		std::cout << "constraint satisfied" << std::endl;
	else
		std::cout << "constraint not satisfied" << std::endl;

	return 0;
}

int expand_terminal_node_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "Expand terminal node test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	ExpandableIRISParams_t& ExpandableIRISParams = ExpandableIRISConic::getDefaultExpandableIRISParams();
	ExpandableIRISParams.IRISParams.n = 2;
	ExpandableIRISConic expandableIRISConic = ExpandableIRISConic(cObs, ExpandableIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	expandableIRISConic.buildNavGraph(qstart, qtarget);
	expandableIRISConic.addConvexSets(qstart);
	expandableIRISConic.expandTerminalNode(3);
	std::vector<int> nodeKeys = expandableIRISConic.gcs.getNodeKeys();
	std::cout << "Results" << std::endl;
	for (int i = 0; i < expandableIRISConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = expandableIRISConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	std::cout << "NavGraph of convex sets" << std::endl;
	expandableIRISConic.gcs.print();
	std::vector<int> edgeNodeKeys = expandableIRISConic.navGraph.getNodeKeys();
	for (int i = 0; i < expandableIRISConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = expandableIRISConic.navGraph.getNode(edgeNodeKeys[i]);
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
	expandableIRISConic.navGraph.print();
	return 0;
}


int addConvexSets_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "Add convex sets test" << std::endl;
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	Eigen::Vector<double, 2> q0({ 0.,0. });
	Eigen::Vector<double, 2> q1({ 4.,-8. });
	Eigen::Vector<double, 2> q2({ 4.5,-5. });
	CObsConic cObs = getCObsPolyhedronV();
	ExpandableIRISParams_t& ExpandableIRISParams = ExpandableIRISConic::getDefaultExpandableIRISParams();
	ExpandableIRISParams.IRISParams.n = 2;
	ExpandableIRISConic expandableIRISConic = ExpandableIRISConic(cObs, ExpandableIRISParams);
	expandableIRISConic.buildNavGraph(qstart, qtarget);
	std::cout << "Adding seed 0 " << std::endl;
	std::cout << q0 << std::endl;
	expandableIRISConic.addConvexSets(q0);
	std::cout << "Adding seed 1 " << std::endl;
	std::cout << q1 << std::endl;
	expandableIRISConic.addConvexSets(q1);
	std::cout << "Adding seed 2 " << std::endl;
	std::cout << q2 << std::endl;
	expandableIRISConic.addConvexSets(q2);
	std::vector<int> nodeKeys = expandableIRISConic.gcs.getNodeKeys();
	std::cout << "Results" << std::endl;
	for (int i = 0; i < expandableIRISConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = expandableIRISConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	std::cout << "NavGraph of convex sets" << std::endl;
	expandableIRISConic.gcs.print();
	std::vector<int> edgeNodeKeys = expandableIRISConic.navGraph.getNodeKeys();
	for (int i = 0; i < expandableIRISConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = expandableIRISConic.navGraph.getNode(edgeNodeKeys[i]);
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
	expandableIRISConic.navGraph.print();
	return 0;
}

int addAllConvexSets_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "Add all convex sets test" << std::endl;
	Range* range = Range::getInstance();
	double lb = -20.;
	double ub = 20.;
	range->setRange(2, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_15_obstacles();
	ExpandableIRISParams_t& ExpandableIRISParams = ExpandableIRISConic::getDefaultExpandableIRISParams();
	ExpandableIRISParams.IRISParams.n = 2;
	ExpandableIRISParams.IRISParams.shrinkFactor = 0.125;
	ExpandableIRISParams.IRISParams.maxTrials = 999;
	std::ofstream fout("ExpandableIRIS_generation_test1_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	cObs.print(fout);
	ExpandableIRISConic expandableIRISConic = ExpandableIRISConic(cObs, ExpandableIRISParams);
	expandableIRISConic.buildNavGraph(qstart, qtarget);
	int trials = 0;
	while (trials < ExpandableIRISParams.IRISParams.maxTrials)
	{
		Eigen::VectorXd seed = expandableIRISConic.generateRandomSeed(trials);
		if (trials >= ExpandableIRISParams.IRISParams.maxTrials)
			break;
		expandableIRISConic.addConvexSets(seed,false);  //When we call this function, we allocate memory of the generated convexSets. Where we should delete them?
	}
	std::vector<int> nodeKeys = expandableIRISConic.gcs.getNodeKeys();
	std::cout << "Results" << std::endl;
	for (int i = 0; i < expandableIRISConic.gcs.numNodes; i++)
	{
		std::cout << "Convex set " << nodeKeys[i] << std::endl;
		Node* node = expandableIRISConic.gcs.getNode(nodeKeys[i]);
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		polyNode->polyhedron.print();
	}
	std::cout << "NavGraph of convex sets" << std::endl;
	expandableIRISConic.gcs.print();
	std::vector<int> edgeNodeKeys = expandableIRISConic.navGraph.getNodeKeys();
	for (int i = 0; i < expandableIRISConic.navGraph.numNodes; i++)
	{
		std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
		Node* node = expandableIRISConic.navGraph.getNode(edgeNodeKeys[i]);
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
	expandableIRISConic.navGraph.print();
	expandableIRISConic.gcs.print(fout, trials);
	expandableIRISConic.navGraph.printGraph(fout, trials);

	fout << "figure; " << std::endl;
	fout << "fGraph" << trials << "=plot(g" << trials << ");" << std::endl;

	fout << "allA=[AObs AGCS" << trials << "]; " << std::endl;
	fout << "allB=[bObs bGCS" << trials << "];" << std::endl;
	fout << "allColors=[colObs colGCS" << trials << "];" << std::endl;
	fout << "fGCS" << trials << "=figure; " << std::endl;
	fout << "plotregion(allA, allB, [], [], allColors);" << std::endl;
	fout << "hold on;" << std::endl;
	fout << "xlabel('X [m.]');" << std::endl;
	fout << "ylabel('Y [m.]');" << std::endl;
	fout << "for i = 1:length(centroidGCS" << trials << ")" << std::endl;
	fout << "  text(centroidGCS" << trials << "{i}(1), centroidGCS" << trials << "{i}(2),nameGCS" << trials << "{i}, 'FontSize', 12); " << std::endl;
	fout << "end" << std::endl;
	fout << "axis equal;" << std::endl;
	fout << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
	return 0;
}

int main(int argc, char** argv)
{
	//terminalNode_randomseed_test();
	//qstart_qtarget_test(Eigen::Vector<double, 2>({ 0.,0. }), Eigen::Vector<double, 2>({9.,8}));
	//addConvexSets_test(Eigen::Vector<double, 2>({ 0.,0. }), Eigen::Vector<double, 2>({ 9.,8 }));
	//expand_terminal_node_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
	addAllConvexSets_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
}