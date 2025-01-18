#include <iostream>
#include "Range.h"
#include "Point.h"
#include "Polyhedron.h"
#include "PolyhedronObstacleCircularRobot.h"
#include "Ellipsoid.h"

int isInsideSeparatingHyperplanes()
{
	std::cout << "Polyhedron isInsideSeparatingHyperplanes Test...";
	Eigen::Matrix<double, 5, 2> A({ {0.5547, 0.83205},{-0.164399, -0.986394},{-0.928477, 0.371391},{-0.707107, 0.707107},{0.393919, -0.919145} });
	Eigen::Vector<double, 5> b({ 4.80058, 2.10681, 4.79855, 4.36149, 2.35798 });
	Polyhedron poly(A, b);
	//Check the first separating hyperplane
	Eigen::MatrixXd A1(1, 2);
	A1 << 0.393919, -0.919145;
	Eigen::VectorXd b1(1);
	b1 << 2.07233;
	bool check1=poly.isInsideSeparatingHyperplanes(A1, b1);
	Eigen::MatrixXd A2(2, 2);
	A2 << 0.393919, -0.919145, -0.948683, 0.316228;
	Eigen::VectorXd b2(2);
	b2 << 2.07233, -5.95701;
	bool check2 = poly.isInsideSeparatingHyperplanes(A2, b2);
	//Now, let's reproduce the condition in which we repeat the calculi, but we need to restart again using the same polyhedron
	Eigen::MatrixXd A0b(0, 2);
	Eigen::VectorXd b0b(0);
	bool check0b = poly.isInsideSeparatingHyperplanes(A0b, b0b);
	Eigen::MatrixXd A1b(1, 2);
	A1b << 0.393919, -0.919145;
	Eigen::VectorXd b1b(1);
	b1b << 2.07233;
	bool check1b = poly.isInsideSeparatingHyperplanes(A1b, b1b);
	Eigen::MatrixXd A2b(2, 2);
	A2b << 0.393919, -0.919145, -0.948683, 0.316228;
	Eigen::VectorXd b2b(2);
	b2b << 2.07233, -6.95701;
	bool check2b = poly.isInsideSeparatingHyperplanes(A2b, b2b);
	
	if (check1 && check2 && check0b && check1b && !check2b)
		std::cout << "Test passed successfully" << std::endl;
	else
		std::cout << "Test failed!" << std::endl;
	return 0;
}

int main(int argc, char** argv)
{
	Range* range = Range::getInstance();
	range->setRange(2, -10., 10.);
	isInsideSeparatingHyperplanes();
	return 0;
}