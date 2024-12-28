#include <iostream>
#include "Range.h"
#include "Point.h"
#include "Polyhedron.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "Ellipsoid.h"

int polyhedronDistanceTest()
{
	double distance_test = 2.0;
	Eigen::Vector<double, 2> p_out_test({ -1.,0.5 });
	std::cout << "Polyhedron Distance Test...";
	Eigen::Matrix<double, 4, 2> A1({ {1.,0.},{0.,1. },{-1.,0},{0.,-1.} });
	Eigen::Vector<double, 4> b1({ 1.,1.,1.,1. });
	Polyhedron poly(A1, b1);
	Eigen::VectorXd p_in(2);
	Eigen::VectorXd p_out(2);
	double distance;
	p_in << -3, 0.5;
	distance = poly.closestPoint(p_in, p_out);
	//std::cout << "distance: " << std::endl;
	//std::cout << distance << std::endl;
	//std::cout << "Point: " << std::endl;
	//std::cout << p_out << std::endl;
	if (((distance - distance_test) < 1e-3) && (((p_out - p_out_test).norm()) < 1e-3))
	{	
		std::cout << "successful" << std::endl;
		return 0;
	}
	else
	{
		std::cout << "failed" << std::endl;
		return -1;
	}
}

int polyhedronVDistanceTest()
{
	double distance_test = 2.0;
	Eigen::Vector<double, 2> p_out_test({ -1.,0.5 });
	std::cout << "PolyhedronV Distance Test...";
	Eigen::Matrix<double, 2, 4> v({ {-1.,1.,1.,-1.},{-1.,-1.,1.,1.} });
	PolyhedronV poly(v);
	Eigen::VectorXd p_in(2);
	Eigen::VectorXd p_out(2);
	double distance;
	p_in << -3, 0.5;
	distance = poly.closestPoint(p_in, p_out);
	//std::cout << "distance: " << std::endl;
	//std::cout << distance << std::endl;
	//std::cout << "Point: " << std::endl;
	//std::cout << p_out << std::endl;
	if (((distance - distance_test) < 1e-3) && (((p_out - p_out_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		return 0;
	}
	else
	{
		std::cout << "failed" << std::endl;
		return -1;
	}
}

int polyhedronObstacleCircleDistanceTest()
{
	double distance_test = 1.7;
	Eigen::Vector<double, 2> p_out_test({ -1.3,0.5 });
	std::cout << "Polyhedron Obstacle Circular Robot Test...";
	Eigen::Matrix<double, 2, 4> v({{-1.,1.,1.,-1.},{-1.,-1.,1.,1.}});
	double r = 0.3;
	PolyhedronObstacleCircularRobot obs(v, r);
	Eigen::VectorXd p_in(2);
	Eigen::VectorXd p_out(2);
	double distance;
	p_in << -3., 0.5;
	distance=obs.closestPoint(p_in, p_out);
	//std::cout << "distance: " << std::endl;
	//std::cout << distance << std::endl;
	//std::cout << "Point: " << std::endl;
	//std::cout << p_out << std::endl;
	if (((distance - distance_test) < 1e-3) && (((p_out - p_out_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		return 0;
	}
	else
	{
		std::cout << "failed" << std::endl;
		return -1;
	}
}

int ellipsoidDistanceTest()
{
	double distance_test = 1.;
	Eigen::Vector<double, 2> p_out_test({ -1.,0.});
	std::cout << "Ellipsoid Test...";
	Eigen::Matrix<double, 2, 2> C1({ {1.,0.},{0., 1.} });
	Eigen::Vector<double, 2> d1({ 0.,0.});
	Ellipsoid ellipsoid1 = Ellipsoid(C1, d1);

	Eigen::VectorXd p_in(2);
	Eigen::VectorXd p_out(2);
	double distance;
	p_in << -2., 0.;
	distance = ellipsoid1.closestPoint(p_in, p_out);
	//std::cout << "Distance to ellipsoid 1 " << std::endl;
	//std::cout << distance << std::endl;
	//std::cout << p_out << std::endl;
	if (((distance - distance_test) < 1e-3) && (((p_out - p_out_test).norm()) < 1e-3))
	{
		std::cout << "successful" << std::endl;
		return 0;
	}
	else
	{
		std::cout << "failed" << std::endl;
		return -1;
	}
}

int boundingBoxTest()
{
	Eigen::Matrix<double, 2, 2> bb0_test({{-2.53852,-1.46148},{-0.538516,0.538516}});
	Eigen::Matrix<double, 2, 2> bb1_test({ {-1.,1.},{-1.,1.}});
	Eigen::Matrix<double, 2, 2> bb2_test({ {-5.,-4.},{-1.,1.} });
	Eigen::Matrix<double, 2, 2> bb3_test({ {-3.3,-0.7},{1.2,2.3}});
	std::cout << "Bounding Box Test..." << std::endl;
	std::cout << "Ellipsoid" << std::endl;
	Eigen::Matrix<double, 2, 2> C({ {0.5,0.2},{0.2,0.5} });
	Eigen::Vector<double, 2> d({ -2.,0. });
	Ellipsoid ellipsoid = Ellipsoid(C, d);
	Eigen::MatrixXd bb0 = ellipsoid.getBoundingBox();
	bool test0 = bb0.isApprox(bb0_test,0.001);
	if (test0)
	{
		std::cout << "successful" << std::endl;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	std::cout << "Polyhedron hyperplanes" << std::endl;
	Eigen::Matrix<double, 4, 2> A1({ {1.,0.},{0.,1.},{-1.,0},{0.,-1.} });
	Eigen::Vector<double, 4> b1({ 1.,1.,1.,1. });
	Polyhedron poly1(A1, b1);
	Eigen::MatrixXd bb1=poly1.getBoundingBox();
	bool test1 = bb1.isApprox(bb1_test, 0.001);
	if (test1)
	{
		std::cout << "successful" << std::endl;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	std::cout << "Polyhedron vertices" << std::endl;
	Eigen::Matrix<double, 2, 5> v2({ {-5., -5., -4., -4., -4.5},{-1., 1., -1., 1.,0.} });
	PolyhedronV poly2(v2);
	Eigen::MatrixXd bb2 = poly2.getBoundingBox();
	bool test2 = bb2.isApprox(bb2_test, 0.001);
	if (test2)
	{
		std::cout << "successful" << std::endl;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	std::cout << "Polyhedron obstacle circular robot" << std::endl;
	Eigen::Matrix<double, 2, 4> v({ {-3.,-1.,-1.,-3.},{1.5,1.5,2.,2.} });
	double r = 0.3;
	PolyhedronObstacleCircularRobot obs(v, r);
	Eigen::MatrixXd bb3 = obs.getBoundingBox();
	bool test3 = bb3.isApprox(bb3_test, 0.001);
	if (test3)
	{
		std::cout << "successful" << std::endl;
	}
	else
	{
		std::cout << "failed" << std::endl;
	}
	return test0 && test1 && test2 && test3;
}

int main(int argc, char** argv)
{
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	polyhedronDistanceTest();
	polyhedronVDistanceTest();
	polyhedronObstacleCircleDistanceTest();
	ellipsoidDistanceTest();
	boundingBoxTest();
	return 0;
}