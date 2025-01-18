#include <iostream>
#include "Range.h"
#include "Point.h"
#include "Sphere.h"
#include "Polyhedron.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "Ellipsoid.h"

int pointExpandingEllipsoidTest()
{
	double distance_test = 5.1287;
	Eigen::Vector<double, 2> p_out_test({ 0.,0. });
	std::cout << "Point Expanding Ellipsoid Distance Test...";
	Eigen::Matrix<double, 2, 2> C1({ {0.5,0.2},{0.2,0.5} });
	Eigen::Vector<double, 2> d1({ -2.,0. });
	Eigen::Vector<double, 2> p({ 0.,0. });
	Point point(p);
	Ellipsoid ellipsoid(C1, d1);
	Eigen::VectorXd p_out(2);
	double distance = point.closestPointExpandingEllipsoid(ellipsoid, p_out);
	std::cout << distance << std::endl;
	std::cout << "Closest point: " << std::endl;
	std::cout << p_out<< std::endl;
	bool test = false;
	if (((distance - distance_test) < 1e-3) && (((p_out - p_out_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		test = true;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	return test ? 0 : -1;
}

int sphereExpandingEllipsoidTest()
{
	double distance_test = 4.2055;
	Eigen::Vector<double, 2> p_out_test({-0.2529,0.1614});
	std::cout << "Sphere Expanding Ellipsoid Distance Test...";
	Eigen::Matrix<double, 2, 2> C1({ {0.5,0.2},{0.2,0.5} });
	Eigen::Vector<double, 2> d1({ -2.,0. });
	Eigen::Vector<double, 2> p({ 0.,0. });
	double r = 0.3;
	Sphere sphere(p,r);
	sphere.allocateClosestPointEllipsoidSolver();
	Ellipsoid ellipsoid(C1, d1);
	Eigen::VectorXd p_out(2);
	double distance = sphere.closestPointExpandingEllipsoid(ellipsoid, p_out);
	std::cout << distance << std::endl;
	std::cout << "Closest point: " << std::endl;
	std::cout << p_out << std::endl;
	bool test = false;
	if (((distance - distance_test) < 1e-3) && (((p_out - p_out_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		test = true;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	return test ? 0 : -1;
}

int ellipsoidExpandingEllipsoidTest()
{
	double distance_test = 2.8263;
	Eigen::Vector<double, 2> p_out_test({-0.5135,0.7886});
	std::cout << "Ellipsoid Expanding Ellipsoid Distance Test...";
	Eigen::Matrix<double, 2, 2> C1({{0.5,0.2},{0.2,0.5} });
	Eigen::Vector<double, 2> d1({-2.,0.});
	Eigen::Matrix<double, 2, 2> C2({{0.5,-0.2},{-0.2,1.} });
	Eigen::Vector<double, 2> d2({0.,0.});
	Ellipsoid ellipsoid1(C1, d1);
	Ellipsoid ellipsoid2(C2,d2);
	ellipsoid2.allocateClosestPointEllipsoidSolver();
	Eigen::VectorXd p_out(2);
	double distance = ellipsoid2.closestPointExpandingEllipsoid(ellipsoid1, p_out);
	std::cout << distance << std::endl;
	std::cout << "Closest point: " << std::endl;
	std::cout << p_out << std::endl;
	bool test = false;
	if (((distance - distance_test) < 1e-3) && (((p_out - p_out_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		test = true;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	return test ? 0 : -1;
}

int polyhedronExpandingEllipsoidTest()
{
	double distance1_test = 1.857;
	Eigen::Vector<double, 2> p_out1_test({-1.,0.6897});
	double distance2_test = 3.8392;
	Eigen::Vector<double, 2> p_out2_test({-4.,-1.});
	double distance3_test = 2.2283;
	Eigen::Vector<double, 2> p_out3_test({-1.1724,1.2});
	std::cout << "Polyhedron Expanding Ellipsoid Distance Test..." << std::endl;
	Eigen::Matrix<double, 4, 2> A1({ {1.,0.},{0.,1.},{-1.,0},{0.,-1.} });
	Eigen::Vector<double, 4> b1({ 1.,1.,1.,1. });
	Polyhedron poly1(A1, b1);
	//poly1.allocateClosestPointEllipsoidSolver();
	Eigen::Matrix<double, 2, 5> v2({ {-5., -5., -4., -4., -4.5},{-1., 1., -1., 1.,0.}});
	PolyhedronV poly2(v2);
	//poly2.allocateClosestPointEllipsoidSolver();
	Eigen::Matrix<double, 2, 4> v({{-3.,-1.,-1.,-3.},{1.5,1.5,2.,2.}});
	double r = 0.3;
	PolyhedronObstacleCircularRobot obs(v, r);
	//obs.allocateClosestPointEllipsoidSolver();
	Eigen::Matrix<double, 2, 2> C1({ {0.5,0.2},{0.2,0.5}});
	Eigen::Vector<double, 2> d1({-2.,0.});
	Ellipsoid ellipsoid1 = Ellipsoid(C1, d1);
	std::cout << "Distance to polyhedron hyperplanes: ";
	double distance1;
	Eigen::VectorXd p_out1(2);
	distance1=poly1.closestPointExpandingEllipsoid(ellipsoid1, p_out1);
	std::cout << distance1 << std::endl;
	std::cout << "Closest point: " << std::endl;
	std::cout << p_out1 << std::endl;
	bool test1 = false;
	if (((distance1 - distance1_test) < 1e-3) && (((p_out1 - p_out1_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		test1 = true;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	double distance2;
	Eigen::VectorXd p_out2(2);
	distance2 = poly2.closestPointExpandingEllipsoid(ellipsoid1, p_out2);
	std::cout << "Distance to polyhedron vertices: ";
	std::cout << distance2 << std::endl;
	std::cout << "Closest point: " << std::endl;
	std::cout << p_out2 << std::endl;
	bool test2 = false;
	if (((distance2 - distance2_test) < 1e-3) && (((p_out2 - p_out2_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		test2 = true;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	double distance3;
	Eigen::VectorXd p_out3(2);
	distance3 = obs.closestPointExpandingEllipsoid(ellipsoid1, p_out3);
	std::cout << "Distance to polyhedron obstacle circular robot: ";
	std::cout << distance3 << std::endl;
	std::cout << "Closest point: " << std::endl;
	std::cout << p_out3 << std::endl;
	bool test3 = false;
	if (((distance3 - distance3_test) < 1e-3) && (((p_out3 - p_out3_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		test3 = true;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	if (test1 && test2 && test3)
		return 0;
	else
		return -1;
}

int main(int argc, char** argv)
{
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	pointExpandingEllipsoidTest();
	sphereExpandingEllipsoidTest();
	ellipsoidExpandingEllipsoidTest();
	polyhedronExpandingEllipsoidTest();
	return 0;
}