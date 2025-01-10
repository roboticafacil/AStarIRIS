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
#include "ConvexRelaxationMinDistanceSolver.h"
#include "MinDistanceSolver.h"

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
	Eigen::Matrix<double, 2, 3> v3({ {0.,1.,2.},{-8.,-6.,-8.} });
	PolyhedronV* o3 = new PolyhedronV(v3);
	cObs.addObject(o3);
	Eigen::Matrix<double, 2, 5> v4({ {5.,5.,6.,7.,7.},{-7.,-5.,-4.,-5.,-7.} });
	PolyhedronV* o4 = new PolyhedronV(v4);
	cObs.addObject(o4);
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

int main(int argc, char** argv)
{
	terminalNode_randomseed_test();
	qstart_qtarget_test(Eigen::Vector<double, 2>({ 0.,0. }), Eigen::Vector<double, 2>({9.,8}));
	addConvexSets_test(Eigen::Vector<double, 2>({ 0.,0. }), Eigen::Vector<double, 2>({ 9.,8 }));
	expand_terminal_node_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
}