#include <iostream>
#include "Range.h"
#include "Point.h"
#include "Polyhedron.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "Ellipsoid.h"
#include <chrono>

int removeConstraintsTest()
{
	std::cout << "Remove Constraint Test..." << std::endl;
	Eigen::Matrix<double, 4, 2> A1({ {-0.707074,  0.707139},{ 0.900728, -0.434384},{  0.39392, -0.919145},{ 0.397836,  0.917457} });
	Eigen::Vector<double, 4> b1({ -1.41438,0.576544, 3.43038,-1.46346});
	Polyhedron poly1(A1, b1);
	std::cout << "Original polyhedron" << std::endl;
	poly1.print();
	Eigen::MatrixXd A2;
	Eigen::VectorXd b2;
	std::chrono::duration<double, std::milli> ms_double;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	t1 = std::chrono::high_resolution_clock::now();
	Polyhedron::removeConstraints(A1, b1, A2, b2);
	Polyhedron poly2(A2, b2);
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	std::cout << "Updated polyhedron" << std::endl;
	poly2.print();
	std::cout << "Static constraint removal performance" << std::endl;
	std::cout << ms_double.count() << "ms" << std::endl;
	return 0;
}

int removeConstraintsInplaceTest()
{
	std::cout << "Remove Constraint (in-place) Test..." << std::endl;
	Eigen::Matrix<double, 4, 2> A1({ {-0.707074,  0.707139},{ 0.900728, -0.434384},{  0.39392, -0.919145},{ 0.397836,  0.917457} });
	Eigen::Vector<double, 4> b1({ -1.41438,0.576544, 3.43038,-1.46346 });
	Polyhedron poly(A1, b1);
	std::cout << "Original polyhedron" << std::endl;
	poly.print();
	std::chrono::duration<double, std::milli> ms_double;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	t1 = std::chrono::high_resolution_clock::now();
	poly.removeConstraints();
	t2 = std::chrono::high_resolution_clock::now();
	ms_double = t2 - t1;
	std::cout << "Updated polyhedron" << std::endl;
	poly.print();
	std::cout << "In-place constraint removal performance" << std::endl;
	std::cout << ms_double.count() << "ms" << std::endl;
	return 0;
}

int main(int argc, char** argv)
{
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	for (int i=0;i<10;i++)
		removeConstraintsTest();
	for (int i = 0; i < 10; i++)
		removeConstraintsInplaceTest();
	return 0;
}