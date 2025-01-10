#include<Eigen/Dense>
#include <chrono>
#include "Range.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "GCS.h"
#include "NavGraph.h"
#include "CObsConic.h"
#include "IRISConic.h"
#include "MCIRISConic.h"
#include "PolyhedronNode.h"
#include "PointNode.h"
#include "ConvexRelaxationMinDistanceSolver.h"
#include <fstream>

CObsConic getCObsPolyhedronV_5_obstacles()
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

CObsConic getCObsPolyhedronV_16_obstacles()
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
	Eigen::Matrix<double, 2, 4> v15({ {0., 8., 8., 0.},{12., 12., 13., 13.} });
	PolyhedronV* o15 = new PolyhedronV(v15);
	cObs.addObject(o14);
	return cObs;
}


int MC_IRIS_test(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "MC IRIS test" << std::endl;
	Range* range = Range::getInstance();
	double lb = -10.;
	double ub = 10.;
	range->setRange(2, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_5_obstacles();
	MCIRISParams_t& MCIRISParams = MCIRISConic::getDefaultMCIRISParams();
	MCIRISParams.ExpandableIRISParams.maxItersOptimalPath = 1;
	MCIRISParams.maxIters = 100;
	MCIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	std::ofstream fout("MC_IRIS_test.m");
	fout << "close all;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	cObs.print(fout);
	
	MCIRISConic irisConic = MCIRISConic(cObs, MCIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.do_MCRelaxedSolver();	
	irisConic.gcs.print(fout);
	fout << "allA=[AObs AGCS];" << std::endl;
	fout << "allB=[bObs bGCS];" << std::endl;
	fout << "allColors=[colObs colGCS];" << std::endl;
	fout << "figure;" << std::endl;
	fout << "plotregion(allA, allB, [], [], allColors);" << std::endl;
	fout << "hold on;" << std::endl;
	fout << "xlabel('X [m.]');" << std::endl;
	fout << "ylabel('Y [m.]');" << std::endl;
	fout << "for i = 1:length(centroidGCS)" << std::endl;
	fout << "  text(centroidGCS{i}(1), centroidGCS{i}(2),nameGCS{i},'FontSize',12);" << std::endl;
	fout << "end" << std::endl;
	fout << "axis equal;" << std::endl;
	fout << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
	
	std::cout << "Feasible Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << irisConic.feasibleSolution.x << std::endl;
	fout << "FeasibleSol=[" << irisConic.feasibleSolution.x << "];" << std::endl;
	fout << "plot(FeasibleSol(:, 1), FeasibleSol(:, 2), 'b', 'LineWidth', 3);" << std::endl;
	fout << "plot(FeasibleSol(:, 1), FeasibleSol(:, 2), 'b.', 'MarkerSize', 12, 'LineWidth', 3);" << std::endl;
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = irisConic.optimalPath.nodeKeys.begin(); it != irisConic.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = irisConic.optimalPath.edgeKeys.begin(); it != irisConic.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << irisConic.feasibleSolution.cost << std::endl;
	std::cout << "Results writen to file successfully" << std::endl; 
	return 0;
}

int MC_IRIS_test1(Eigen::Vector<double, 2>& qstart, Eigen::Vector<double, 2>& qtarget)
{
	std::cout << "MC IRIS test1" << std::endl;
	Range* range = Range::getInstance();
	double lb = -20.;
	double ub = 20.;
	range->setRange(2, lb, ub);
	CObsConic cObs = getCObsPolyhedronV_16_obstacles();
	std::ofstream fout("MC_IRIS_test1.m");
	fout << "close all;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;
	cObs.print(fout);
	MCIRISParams_t& MCIRISParams = MCIRISConic::getDefaultMCIRISParams();
	MCIRISParams.ExpandableIRISParams.maxItersOptimalPath = 1;
	MCIRISParams.maxIters = 25;
	MCIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	MCIRISConic irisConic = MCIRISConic(cObs, MCIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.do_MCRelaxedSolver();
	irisConic.gcs.print(fout);
	fout << "allA=[AObs AGCS];" << std::endl;
	fout << "allB=[bObs bGCS];" << std::endl;
	fout << "allColors=[colObs colGCS];" << std::endl;
	fout << "figure;" << std::endl;
	fout << "plotregion(allA, allB, [], [], allColors);" << std::endl;
	fout << "hold on;" << std::endl;
	fout << "xlabel('X [m.]');" << std::endl;
	fout << "ylabel('Y [m.]');" << std::endl;
	fout << "for i = 1:length(centroidGCS)" << std::endl;
	fout << "  text(centroidGCS{i}(1), centroidGCS{i}(2),nameGCS{i},'FontSize',12);" << std::endl;
	fout << "end" << std::endl;
	fout << "axis equal;" << std::endl;
	fout << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;

	std::cout << "Feasible Solution " << std::endl;
	std::cout << "x variables " << std::endl;
	std::cout << irisConic.feasibleSolution.x << std::endl;
	fout << "FeasibleSol=[" << irisConic.feasibleSolution.x << "];" << std::endl;
	fout << "plot(FeasibleSol(:, 1), FeasibleSol(:, 2), 'b', 'LineWidth', 3);" << std::endl;
	fout << "plot(FeasibleSol(:, 1), FeasibleSol(:, 2), 'b.', 'MarkerSize', 12, 'LineWidth', 3);" << std::endl;
	std::cout << "nodes keys ";
	for (std::vector<int>::iterator it = irisConic.optimalPath.nodeKeys.begin(); it != irisConic.optimalPath.nodeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "edge keys ";
	for (std::vector<int>::iterator it = irisConic.optimalPath.edgeKeys.begin(); it != irisConic.optimalPath.edgeKeys.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
	std::cout << "cost " << std::endl;
	std::cout << irisConic.feasibleSolution.cost << std::endl;
	std::cout << "Results writen to file successfully" << std::endl;
	return 0;
}

int main(int argc, char** argv)
{
	//MC_IRIS_test(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
	MC_IRIS_test1(Eigen::Vector<double, 2>({ -5.,-8. }), Eigen::Vector<double, 2>({ 9.,8 }));
}