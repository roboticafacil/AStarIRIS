#include<Eigen/Dense>
#include <chrono>
#include "Range.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "GCS.h"
#include "CObsConic.h"
#include "IRISConic.h"
#include "ExpandableIRISConic.h"
#include "PolyhedronNode.h"
#include "AStarIRISConic.h"
#include <fstream>

CObsConic getCObsPolyhedronV_1_obstacles4D()
{
	CObsConic cObs = CObsConic();
	Eigen::Matrix<double, 4, 4> A0({ {0.000000,1.000000,0.000000,0.000000},{-0.000000,-1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-1.000000,-0.000000,0.000000,0.000000} });
	Eigen::Vector<double, 4> b0({ 3.000000,3.000000,0.000000,10.000000 });
	Polyhedron* o0 = new Polyhedron(A0, b0);
	cObs.addObject(o0);
	Eigen::Matrix<double, 4, 4> A1({ {0.000000,0.000000,0.000000,1.000000},{0.000000,0.000000,-0.000000,-1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-1.000000,-0.000000} });
	Eigen::Vector<double, 4> b1({ 3.000000,3.000000,0.000000,10.000000 });
	Polyhedron* o1 = new Polyhedron(A1, b1);
	cObs.addObject(o1);
	return cObs;
}

CObsConic getCObsPolyhedronV_15_obstacles4D()
{
	CObsConic cObs = CObsConic();
	Eigen::Matrix<double, 4, 4> A0({ {0.707107,0.707107,0.000000,0.000000},{0.707107,-0.707107,0.000000,0.000000},{-0.707107,0.707107,0.000000,0.000000},{-0.707107,-0.707107,0.000000,0.000000} });
	Eigen::Vector<double, 4> b0({ 1.414214,-5.656854,7.071068,0.000000 });
	Polyhedron* o0 = new Polyhedron(A0, b0);
	cObs.addObject(o0);
	Eigen::Matrix<double, 4, 4> A1({ {0.000000,0.000000,0.707107,0.707107},{0.000000,0.000000,0.707107,-0.707107},{0.000000,0.000000,-0.707107,0.707107},{0.000000,0.000000,-0.707107,-0.707107} });
	Eigen::Vector<double, 4> b1({ 1.414214,-5.656854,7.071068,0.000000 });
	Polyhedron* o1 = new Polyhedron(A1, b1);
	cObs.addObject(o1);
	Eigen::Matrix<double, 5, 4> A2({ {-0.554700,-0.832050,0.000000,0.000000},{0.447214,0.894427,0.000000,0.000000},{0.948683,-0.316228,0.000000,0.000000},{0.894427,0.447214,0.000000,0.000000},{-0.970143,0.242536,0.000000,0.000000} });
	Eigen::Vector<double, 5> b2({ -6.101702,9.838699,6.957011,10.285913,-3.880570 });
	Polyhedron* o2 = new Polyhedron(A2, b2);
	cObs.addObject(o2);
	Eigen::Matrix<double, 5, 4> A3({ {0.000000,0.000000,-0.554700,-0.832050},{0.000000,0.000000,0.447214,0.894427},{0.000000,0.000000,0.948683,-0.316228},{0.000000,0.000000,0.894427,0.447214},{0.000000,0.000000,-0.970143,0.242536} });
	Eigen::Vector<double, 5> b3({ -6.101702,9.838699,6.957011,10.285913,-3.880570 });
	Polyhedron* o3 = new Polyhedron(A3, b3);
	cObs.addObject(o3);
	Eigen::Matrix<double, 4, 4> A4({ {0.164399,0.986394,0.000000,0.000000},{0.707107,-0.707107,0.000000,0.000000},{0.393919,-0.919145,0.000000,0.000000},{-0.928477,0.371391,0.000000,0.000000} });
	Eigen::Vector<double, 4> b4({ -3.123581,1.414214,2.888742,5.756555 });
	Polyhedron* o4 = new Polyhedron(A4, b4);
	cObs.addObject(o4);
	Eigen::Matrix<double, 4, 4> A5({ {0.000000,0.000000,0.164399,0.986394},{0.000000,0.000000,0.707107,-0.707107},{0.000000,0.000000,0.393919,-0.919145},{0.000000,0.000000,-0.928477,0.371391} });
	Eigen::Vector<double, 4> b5({ -3.123581,1.414214,2.888742,5.756555 });
	Polyhedron* o5 = new Polyhedron(A5, b5);
	cObs.addObject(o5);
	Eigen::Matrix<double, 3, 4> A6({ {-0.000000,-1.000000,0.000000,0.000000},{0.894427,0.447214,0.000000,0.000000},{-0.894427,0.447214,0.000000,0.000000} });
	Eigen::Vector<double, 3> b6({ 9.000000,-2.236068,-4.024922 });
	Polyhedron* o6 = new Polyhedron(A6, b6);
	cObs.addObject(o6);
	Eigen::Matrix<double, 3, 4> A7({ {0.000000,0.000000,-0.000000,-1.000000},{0.000000,0.000000,0.894427,0.447214},{0.000000,0.000000,-0.894427,0.447214} });
	Eigen::Vector<double, 3> b7({ 9.000000,-2.236068,-4.024922 });
	Polyhedron* o7 = new Polyhedron(A7, b7);
	cObs.addObject(o7);
	Eigen::Matrix<double, 5, 4> A8({ {-0.000000,-1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-1.000000,-0.000000,0.000000,0.000000},{0.707107,0.707107,0.000000,0.000000},{-0.707107,0.707107,0.000000,0.000000} });
	Eigen::Vector<double, 5> b8({ 7.000000,7.000000,-5.000000,1.414214,-7.071068 });
	Polyhedron* o8 = new Polyhedron(A8, b8);
	cObs.addObject(o8);
	Eigen::Matrix<double, 5, 4> A9({ {0.000000,0.000000,-0.000000,-1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-1.000000,-0.000000},{0.000000,0.000000,0.707107,0.707107},{0.000000,0.000000,-0.707107,0.707107} });
	Eigen::Vector<double, 5> b9({ 7.000000,7.000000,-5.000000,1.414214,-7.071068 });
	Polyhedron* o9 = new Polyhedron(A9, b9);
	cObs.addObject(o9);
	Eigen::Matrix<double, 5, 4> A10({ {0.948683,0.316228,0.000000,0.000000},{0.658505,0.752577,0.000000,0.000000},{-0.792624,-0.609711,0.000000,0.000000},{0.707107,-0.707107,0.000000,0.000000},{-0.707107,0.707107,0.000000,0.000000} });
	Eigen::Vector<double, 5> b10({ -12.649111,-14.110813,18.901034,5.656854,10.606602 });
	Polyhedron* o10 = new Polyhedron(A10, b10);
	cObs.addObject(o10);
	Eigen::Matrix<double, 5, 4> A11({ {0.000000,0.000000,0.948683,0.316228},{0.000000,0.000000,0.658505,0.752577},{0.000000,0.000000,-0.792624,-0.609711},{0.000000,0.000000,0.707107,-0.707107},{0.000000,0.000000,-0.707107,0.707107} });
	Eigen::Vector<double, 5> b11({ -12.649111,-14.110813,18.901034,5.656854,10.606602 });
	Polyhedron* o11 = new Polyhedron(A11, b11);
	cObs.addObject(o11);
	Eigen::Matrix<double, 4, 4> A12({ {0.000000,1.000000,0.000000,0.000000},{-0.000000,-1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-1.000000,-0.000000,0.000000,0.000000} });
	Eigen::Vector<double, 4> b12({ -13.000000,16.000000,15.000000,0.000000 });
	Polyhedron* o12 = new Polyhedron(A12, b12);
	cObs.addObject(o12);
	Eigen::Matrix<double, 4, 4> A13({ {0.000000,0.000000,0.000000,1.000000},{0.000000,0.000000,-0.000000,-1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-1.000000,-0.000000} });
	Eigen::Vector<double, 4> b13({ -13.000000,16.000000,15.000000,0.000000 });
	Polyhedron* o13 = new Polyhedron(A13, b13);
	cObs.addObject(o13);
	Eigen::Matrix<double, 4, 4> A14({ {-0.000000,1.000000,0.000000,0.000000},{0.000000,-1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-0.707107,-0.707107,0.000000,0.000000} });
	Eigen::Vector<double, 4> b14({ -10.000000,11.000000,-5.000000,12.020815 });
	Polyhedron* o14 = new Polyhedron(A14, b14);
	cObs.addObject(o14);
	Eigen::Matrix<double, 4, 4> A15({ {0.000000,0.000000,-0.000000,1.000000},{0.000000,0.000000,0.000000,-1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-0.707107,-0.707107} });
	Eigen::Vector<double, 4> b15({ -10.000000,11.000000,-5.000000,12.020815 });
	Polyhedron* o15 = new Polyhedron(A15, b15);
	cObs.addObject(o15);
	Eigen::Matrix<double, 4, 4> A16({ {0.000000,1.000000,0.000000,0.000000},{-0.000000,-1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-1.000000,-0.000000,0.000000,0.000000} });
	Eigen::Vector<double, 4> b16({ -5.000000,10.000000,16.000000,-15.000000 });
	Polyhedron* o16 = new Polyhedron(A16, b16);
	cObs.addObject(o16);
	Eigen::Matrix<double, 4, 4> A17({ {0.000000,0.000000,0.000000,1.000000},{0.000000,0.000000,-0.000000,-1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-1.000000,-0.000000} });
	Eigen::Vector<double, 4> b17({ -5.000000,10.000000,16.000000,-15.000000 });
	Polyhedron* o17 = new Polyhedron(A17, b17);
	cObs.addObject(o17);
	Eigen::Matrix<double, 4, 4> A18({ {0.000000,1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-1.000000,-0.000000,0.000000,0.000000},{0.894427,-0.447214,0.000000,0.000000} });
	Eigen::Vector<double, 4> b18({ 15.000000,13.000000,-12.000000,9.391486 });
	Polyhedron* o18 = new Polyhedron(A18, b18);
	cObs.addObject(o18);
	Eigen::Matrix<double, 4, 4> A19({ {0.000000,0.000000,0.000000,1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-1.000000,-0.000000},{0.000000,0.000000,0.894427,-0.447214} });
	Eigen::Vector<double, 4> b19({ 15.000000,13.000000,-12.000000,9.391486 });
	Polyhedron* o19 = new Polyhedron(A19, b19);
	cObs.addObject(o19);
	Eigen::Matrix<double, 4, 4> A20({ {0.000000,1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-1.000000,-0.000000,0.000000,0.000000},{-0.894427,-0.447214,0.000000,0.000000} });
	Eigen::Vector<double, 4> b20({ 15.000000,-12.000000,13.000000,9.391486 });
	Polyhedron* o20 = new Polyhedron(A20, b20);
	cObs.addObject(o20);
	Eigen::Matrix<double, 4, 4> A21({ {0.000000,0.000000,0.000000,1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-1.000000,-0.000000},{0.000000,0.000000,-0.894427,-0.447214} });
	Eigen::Vector<double, 4> b21({ 15.000000,-12.000000,13.000000,9.391486 });
	Polyhedron* o21 = new Polyhedron(A21, b21);
	cObs.addObject(o21);
	Eigen::Matrix<double, 4, 4> A22({ {0.000000,1.000000,0.000000,0.000000},{-0.000000,-1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-1.000000,-0.000000,0.000000,0.000000} });
	Eigen::Vector<double, 4> b22({ 10.000000,-9.000000,5.000000,5.000000 });
	Polyhedron* o22 = new Polyhedron(A22, b22);
	cObs.addObject(o22);
	Eigen::Matrix<double, 4, 4> A23({ {0.000000,0.000000,0.000000,1.000000},{0.000000,0.000000,-0.000000,-1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-1.000000,-0.000000} });
	Eigen::Vector<double, 4> b23({ 10.000000,-9.000000,5.000000,5.000000 });
	Polyhedron* o23 = new Polyhedron(A23, b23);
	cObs.addObject(o23);
	Eigen::Matrix<double, 4, 4> A24({ {0.000000,1.000000,0.000000,0.000000},{-0.000000,-1.000000,0.000000,0.000000},{1.000000,0.000000,0.000000,0.000000},{-1.000000,-0.000000,0.000000,0.000000} });
	Eigen::Vector<double, 4> b24({ 2.000000,5.000000,10.000000,-9.000000 });
	Polyhedron* o24 = new Polyhedron(A24, b24);
	cObs.addObject(o24);
	Eigen::Matrix<double, 4, 4> A25({ {0.000000,0.000000,0.000000,1.000000},{0.000000,0.000000,-0.000000,-1.000000},{0.000000,0.000000,1.000000,0.000000},{0.000000,0.000000,-1.000000,-0.000000} });
	Eigen::Vector<double, 4> b25({ 2.000000,5.000000,10.000000,-9.000000 });
	Polyhedron* o25 = new Polyhedron(A25, b25);
	cObs.addObject(o25);
	Eigen::Matrix<double, 4, 4> A26({ {0.707107,0.707107,0.000000,0.000000},{0.707107,-0.707107,0.000000,0.000000},{-0.707107,0.707107,0.000000,0.000000},{-0.707107,-0.707107,0.000000,0.000000} });
	Eigen::Vector<double, 4> b26({ -0.707107,-9.192388,10.606602,2.121320 });
	Polyhedron* o26 = new Polyhedron(A26, b26);
	cObs.addObject(o26);
	Eigen::Matrix<double, 4, 4> A27({ {0.000000,0.000000,0.707107,0.707107},{0.000000,0.000000,0.707107,-0.707107},{0.000000,0.000000,-0.707107,0.707107},{0.000000,0.000000,-0.707107,-0.707107} });
	Eigen::Vector<double, 4> b27({ -0.707107,-9.192388,10.606602,2.121320 });
	Polyhedron* o27 = new Polyhedron(A27, b27);
	cObs.addObject(o27);
	Eigen::Matrix<double, 4, 4> A28({ {0.707107,0.707107,0.000000,0.000000},{0.707107,-0.707107,0.000000,0.000000},{-0.707107,0.707107,0.000000,0.000000},{-0.707107,-0.707107,0.000000,0.000000} });
	Eigen::Vector<double, 4> b28({ -9.192388,-7.778175,9.192388,10.606602 });
	Polyhedron* o28 = new Polyhedron(A28, b28);
	cObs.addObject(o28);
	Eigen::Matrix<double, 4, 4> A29({ {0.000000,0.000000,0.707107,0.707107},{0.000000,0.000000,0.707107,-0.707107},{0.000000,0.000000,-0.707107,0.707107},{0.000000,0.000000,-0.707107,-0.707107} });
	Eigen::Vector<double, 4> b29({ -9.192388,-7.778175,9.192388,10.606602 });
	Polyhedron* o29 = new Polyhedron(A29, b29);
	cObs.addObject(o29);
	return cObs;
}

int IRIS_test(const double& shrinkFactor, const int& maxTrials, const double& tol)
{
	std::cout << "Multiagent IRIS test" << std::endl;
	Range* range = Range::getInstance();
	double lb = -10.;
	double ub = 10.;
	int agents = 2;
	range->setRange(4, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_1_obstacles4D();
	IRISParams_t& IRISParams = IRISConic::getDefaultIRISParams();
	IRISParams.n = 4;
	IRISParams.shrinkFactor = shrinkFactor;
	IRISParams.maxTrials = maxTrials;
	IRISParams.tol = tol;
	IRISParams.seperatingHyperplaneAligned = true;
	std::ofstream fout("Multiagent_IRIS_generation_test_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	//cObs.print(fout);
	fout << "obstacles1;" << std::endl;
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
	irisConic.gcs.print(fout, 1);
	if (irisConic.gcs.numEdges > 0)
	{
		irisConic.gcs.printGraph(fout);

		fout << "figure; " << std::endl;
		fout << "fGraph1=plot(g1);" << std::endl;
	}
	irisConic.gcs.print(fout,1,agents);
	for (int agent = 1; agent <= agents; agent++)
	{
		fout << "allA=[AObs AGCS1{:," << agent << "}]; " << std::endl;
		fout << "allB=[bObs bGCS1{:," << agent << "}];" << std::endl;
		fout << "allColors=[colObs colGCS1]; " << std::endl;
		fout << "fGCS1=figure; " << std::endl;
		fout << "plotregion(allA, allB, [], [], allColors);" << std::endl;
		fout << "hold on;" << std::endl;
		fout << "xlabel('X [m.]');" << std::endl;
		fout << "ylabel('Y [m.]');" << std::endl;
		fout << "for i = 1:size(centroidGCS1,1)" << std::endl;
		fout << "  text(centroidGCS1{i," << agent << "}(1), centroidGCS1{i," << agent << "}(2), nameGCS1{i}, 'FontSize', 12); " << std::endl;
		fout << "end" << std::endl;
		fout << "axis equal;" << std::endl;
		fout << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
	}
	return 0;
}


int IRIS_test1(const double& shrinkFactor, const int& maxTrials, const double& tol)
{
	std::cout << "Multiagent IRIS test1" << std::endl;
	Range* range = Range::getInstance();
	double lb = -20.;
	double ub = 20.;
	range->setRange(4, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_15_obstacles4D();
	IRISParams_t& IRISParams = IRISConic::getDefaultIRISParams();
	IRISParams.n = 4;
	IRISParams.shrinkFactor = shrinkFactor;
	IRISParams.maxTrials = maxTrials;
	IRISParams.tol = tol;
	IRISParams.seperatingHyperplaneAligned = true;
	std::ofstream fout("Multiagent_IRIS_generation_test1_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	fout << "obstacles15;" << std::endl;
	//cObs.print(fout);
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
	irisConic.gcs.print(fout,1);
	if (irisConic.gcs.numEdges > 0)
	{
		irisConic.gcs.printGraph(fout);

		fout << "figure; " << std::endl;
		fout << "fGraph1=plot(g1);" << std::endl;
	}

	fout << "allA=[AObs AGCS1]; " << std::endl;
	fout << "allB=[bObs bGCS1];" << std::endl;
	fout << "allColors=[colObs colGCS1];" << std::endl;
	fout << "fGCS1=figure; " << std::endl;
	fout << "plotregion(allA, allB, [], [], allColors);" << std::endl;
	fout << "hold on;" << std::endl;
	fout << "xlabel('X [m.]');" << std::endl;
	fout << "ylabel('Y [m.]');" << std::endl;
	fout << "for i = 1:length(centroidGCS1)" << std::endl;
	fout << "  text(centroidGCS{i}(1), centroidGCS1{i}(2),nameGCS1{i}, 'FontSize', 12); " << std::endl;
	fout << "end" << std::endl;
	fout << "axis equal;" << std::endl;
	fout << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
	return 0;
}

int Expandable_IRIS_test1(Eigen::Vector<double, 4>& qstart, Eigen::Vector<double, 4>& qtarget)
{
	std::cout << "Generate GCS PolyhedronV test" << std::endl;
	Range* range = Range::getInstance();
	double lb = -20.;
	double ub = 20.;
	ExpandableIRISParams_t& ExpandableIRISParams = ExpandableIRISConic::getDefaultExpandableIRISParams();
	ExpandableIRISParams.IRISParams.n = 4;
	ExpandableIRISParams.IRISParams.shrinkFactor = 0.125;
	ExpandableIRISParams.IRISParams.maxTrials = 999;
	range->setRange(ExpandableIRISParams.IRISParams.n, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_15_obstacles4D();
	std::ofstream fout("Multiagent_ExpandableIRIS_generation_test1_results.m");
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
		expandableIRISConic.addConvexSets(seed, false);
		std::cout << "GCS:: Num. nodes " << expandableIRISConic.gcs.numNodes << " num edges " << expandableIRISConic.gcs.numEdges << std::endl;
		std::cout << "navGraph:: Num. nodes " << expandableIRISConic.navGraph.numNodes << " num edges " << expandableIRISConic.navGraph.numEdges << std::endl;
	}

	std::cout << "Done!!" << std::endl;
	//expandableIRISConic.gcs.print(fout, trials);
	expandableIRISConic.navGraph.printGraph(fout, trials);

	fout << "figure; " << std::endl;
	fout << "fGraph" << trials << "=plot(g" << trials << ");" << std::endl;

	/*fout << "allA=[AObs AGCS" << trials << "]; " << std::endl;
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
	fout << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;*/
	return 0;
}

int AStar_IRIS_relaxed_solver_test(Eigen::Vector<double, 4>& qstart, Eigen::Vector<double, 4>& qtarget)
{
	std::cout << "Multiagent AStar Relaxed Solver test" << std::endl;
	Range* range = Range::getInstance();
	double lb = -10.;
	double ub = 10.;
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 4;
	AStarIRISParams.ExpandableIRISParams.IRISParams.seperatingHyperplaneAligned = true;
	AStarIRISParams.ExpandableIRISParams.maxItersOptimalPath = 200;
	AStartIRISDebugLevel_t debug_level;
	debug_level.gcsGraph = true;
	debug_level.gcsConvexSets = true;
	debug_level.navGraph = false;
	debug_level.navGraphConvexSets = false;
	debug_level.video_frames = false;
	debug_level.phase1 = true;
	debug_level.phase2 = true;
	debug_level.multiagent = 2;
	range->setRange(AStarIRISParams.ExpandableIRISParams.IRISParams.n, lb, ub);
	CObsConic cObs4D = getCObsPolyhedronV_1_obstacles4D();
	AStarIRISConic irisConic = AStarIRISConic(cObs4D, AStarIRISParams);
	std::ofstream fout("Multiagent_AStar_IRIS_relaxed_solver_test_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	fout << "obstacles1;" << std::endl;
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::duration<double, std::milli> ms_double;
	double phase1Time = 0.0;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.buildNavGraph(qstart, qtarget);
	ConvexRelaxationMinDistanceSPP_GCS solver(&irisConic.navGraph, irisConic.qStartNodeNavGraphKey, irisConic.qTargetNodeNavGraphKey);
	irisConic.do_RelaxedSolver(solver,fout, debug_level);
	irisConic.gcs.printGraph(fout, 0);
	fout << "figure; " << std::endl;
	fout << "fGraph0=plot(g0,'EdgeLabel',str2num(num2str(g0.Edges.Weights,3)));" << std::endl;
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	phase1Time = ms_double.count();
	std::cout << "Performance" << std::endl;
	std::cout << phase1Time << "ms" << std::endl;
	std::cout << "Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << irisConic.feasibleSolution.x << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << irisConic.feasibleSolution.cost << std::endl;
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = irisConic.optimalPath.nodeKeys.begin(); it != irisConic.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = irisConic.optimalPath.edgeKeys.begin(); it != irisConic.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	return 0;
}

int AStar_IRIS_relaxed_solver_test1(Eigen::Vector<double, 4>& qstart, Eigen::Vector<double, 4>& qtarget)
{
	std::cout << "Multiagent AStar Relaxed Solver test1" << std::endl;
	Range* range = Range::getInstance();
	double lb = -20.;
	double ub = 20.;
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 4;
	AStarIRISParams.ExpandableIRISParams.IRISParams.seperatingHyperplaneAligned = true;
	AStarIRISParams.ExpandableIRISParams.maxItersOptimalPath = 200;
	AStartIRISDebugLevel_t debug_level;
	debug_level.gcsGraph = true;
	debug_level.gcsConvexSets = true;
	debug_level.navGraph = false;
	debug_level.navGraphConvexSets = false;
	debug_level.video_frames = false;
	debug_level.phase1 = true;
	debug_level.phase2 = true;
	debug_level.multiagent = 2;
	range->setRange(AStarIRISParams.ExpandableIRISParams.IRISParams.n, lb, ub);
	CObsConic cObs4D = getCObsPolyhedronV_15_obstacles4D();
	AStarIRISConic irisConic = AStarIRISConic(cObs4D, AStarIRISParams);
	//std::ofstream fout("AStar_IRIS_relaxed_solver_test1_results.m");
	std::ofstream fout("Multiagent_AStar_IRIS_relaxed_solver_test1_results.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	//cObs2D.print(fout);
	fout << "obstacles15;" << std::endl;
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::duration<double, std::milli> ms_double;
	double phase1Time = 0.0;
	t1 = std::chrono::high_resolution_clock::now();
	irisConic.buildNavGraph(qstart, qtarget);
	ConvexRelaxationMinDistanceSPP_GCS solver(&irisConic.navGraph, irisConic.qStartNodeNavGraphKey, irisConic.qTargetNodeNavGraphKey);
	irisConic.do_RelaxedSolver(solver,fout,debug_level);
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	phase1Time = ms_double.count();
	std::cout << "Performance" << std::endl;
	std::cout << phase1Time << "ms" << std::endl;
	std::cout << "Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << irisConic.feasibleSolution.x << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << irisConic.feasibleSolution.cost << std::endl;
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = irisConic.optimalPath.nodeKeys.begin(); it != irisConic.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = irisConic.optimalPath.edgeKeys.begin(); it != irisConic.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	return 0;
}

int main(int argc, char** argv)
{
	//IRIS_test(0.125, 999, 1e-3);
	//IRIS_test1(0.125, 999, 1e-3);
	//Expandable_IRIS_test1(Eigen::Vector<double, 4>({ -5.,-8.,-3.75,-6.75}), Eigen::Vector<double, 4>({ 9.,8.,9.5,6.}));
	AStar_IRIS_relaxed_solver_test(Eigen::Vector<double, 4>({ -5.,-7.5,-4,-8.}), Eigen::Vector<double, 4>({ -5.,7.5,-4.,8.}));
	//AStar_IRIS_relaxed_solver_test1(Eigen::Vector<double, 4>({ -5.,-7.,-4,-7 }), Eigen::Vector<double, 4>({ 9.,8.,9.5,6. }));
}