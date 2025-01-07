#include <cmath>
#include <matplot/matplot.h>
#include "CObsConic.h"
#include<Eigen/Dense>
#include "Range.h"
#include "PolyhedronV.h"
//#include "PolyhedronObstacleCircularRobot.h"
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


using namespace matplot;

CObsConic getCObsPolyhedronV()
{
	CObsConic cObs = CObsConic();
	Eigen::Matrix<double, 2, 4> v0({ {-3.,-4.,-5.,-4.},{5.,6.,5.,4.} });
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

int dump_to_matlab_test()
{
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	Eigen::Matrix<double, 4, 2> A1({ {1.,0.},{0.,1. },{-1.,0},{0.,-1.} });
	Eigen::Vector<double, 4> b1({ 1.,1.,1.,1. });
	Polyhedron poly1(A1, b1);
	poly1.print(std::cout,"A","b",true);
	Eigen::Matrix<double, 2, 4> v0({ {-3.,-4.,-5.,-4.},{5.,6.,5.,4.} });
	PolyhedronV* o0 = new PolyhedronV(v0);
	o0->print(std::cout, "A1", "b1");
	return 0;
}



int fill_test()
{
	CObsConic cObs = getCObsPolyhedronV();
	vector_1d vx;
	vector_1d vy;
	matplot::hold(true);
	for (int i = 0; i < cObs.numObjects; i++)
	{
		PolyhedronV* poly = (PolyhedronV*)cObs.conicObjects[i];
		poly->getFilled2DPolyhedron(vx, vy);
		matplot::fill(vx, vy, "r");
	}
	Eigen::Matrix<double, 4, 2> A1({ {1.,0.},{0.,1. },{-1.,0},{0.,-1.} });
	Eigen::Vector<double, 4> b1({ 1.,1.,1.,1. });
	Polyhedron poly1(A1, b1);
	poly1.getFilled2DPolyhedron(vx, vy);
	matplot::fill(vx, vy, "g");
	//vector_1d vvx({ -1.,1.,1.,-1. });
	//vector_1d vvy({ -1.,-1.,1.,1. });
	//matplot::fill(vvx,vvy, "g");
    //matplot::fill(x, y, red_color);
    axis(equal);
	axis({ -10.,10.,-10.,10. });

    show();
	return 0;
}

int drawConvexSets_test()
{
	Eigen::Vector<double, 2> qstart({ -5.,-8. });
	Eigen::Vector<double, 2> qtarget({ 9.,8 });

	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	CObsConic cObs = getCObsPolyhedronV();
	vector_1d vx;
	vector_1d vy;
	matplot::hold(true);
	for (int i = 0; i < cObs.numObjects; i++)
	{
		PolyhedronV* poly = (PolyhedronV*)cObs.conicObjects[i];
		poly->getFilled2DPolyhedron(vx, vy);
		matplot::fill(vx, vy, "r");
	}
	
	AStarIRISParams_t& AStarIRISParams = AStarIRISConic::getDefaultAStarIRISParams();
	AStarIRISParams.ExpandableIRISParams.IRISParams.n = 2;
	AStarIRISConic irisConic = AStarIRISConic(cObs, AStarIRISParams);
	std::cout << "Creating a convex set from qstart " << std::endl;
	std::cout << qstart << std::endl;
	irisConic.buildNavGraph(qstart, qtarget);
	irisConic.doPhase1_RelaxedSolver();
	std::vector<Node*> nodes=irisConic.gcs.getNodes();
	for (std::vector<Node*>::iterator it=nodes.begin();it!=nodes.end();it++)
	{
		PolyhedronNode* polyNode = (PolyhedronNode*)(*it)->getNodeData();
		polyNode->polyhedron.getFilled2DPolyhedron(vx, vy);
		matplot::fill(vx, vy, "g");
		matplot::plot(vx, vy, "k");
	}
	axis(equal);
	axis({ -10.,10.,-10.,10. });

	show();
	return 0;
}


int main(int argc, char** argv) {
    
	//fill_test();
	dump_to_matlab_test();
	//drawConvexSets_test();
    return 0;
}