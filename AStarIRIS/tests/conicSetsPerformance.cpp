#include <chrono>
#include<Eigen/Dense>
#include <iostream>
#include "Range.h"
#include "Point.h"
#include "Sphere.h"
#include "Ellipsoid.h"
#include "Polyhedron.h"
#include "PolyhedronV.h"
#include "PolyhedronObstacleCircularRobot.h"

void closestPointPerformance()
{
    /*Performance tests for closest point*/
    std::cout << "Closest Point Performance" << std::endl;
    int n = 2;
    int N = 100;
    std::chrono::duration<double, std::milli> ms_double;
    Eigen::MatrixXd points = Eigen::MatrixXd::Random(n, N);
    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;
    double pointTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Point p(points.col(i));
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            p.closestPoint(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            pointTime += ms_double.count();
        }
    }
    std::cout << "Point performance" << std::endl;
    std::cout << pointTime/N << "ms" << std::endl;

    double radius = 0.1;
    double sphereTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Sphere c(points.col(i), radius);
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            c.closestPoint(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            sphereTime += ms_double.count();
        }
    }
    std::cout << "Sphere performance" << std::endl;
    std::cout << sphereTime /N << "ms" << std::endl;

    Eigen::Matrix<double, 2, 2> C({ {0.1,0.035},{0.035,0.2} });
    Eigen::Vector<double, 2> d({ 0.,0. });
    double ellipsoidTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Ellipsoid e(C, points.col(i));
        e.allocateSolver();
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            e.closestPoint(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            ellipsoidTime += ms_double.count();
        }
    }
    std::cout << "Ellipsoid performance" << std::endl;
    std::cout << ellipsoidTime/N << "ms" << std::endl;

    Eigen::Matrix<double, 2, 4> V0({ { -0.1, 0.1, 0.1, -0.1}, {-0.1,-0.1,0.1,0.1} });
    double polyVTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::MatrixXd v = V0.colwise() + points.col(i);
        PolyhedronV poly(v);
        Eigen::VectorXd p_out;
        poly.closestVertex(points.col(i+1), p_out); //Compute once just to ensure the solver is allocated
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            poly.closestVertex(points.col(j), p_out);
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            polyVTime += ms_double.count();
        }
    }
    std::cout << "PolyhedronV performance" << std::endl;
    std::cout << polyVTime/N << "ms" << std::endl;

    Eigen::Matrix<double, 4, 2> A0({ {1.,0},{0.,1.},{-1.,0},{0.,-1.} });
    Eigen::Vector<double, 4> b0({ 0.1,0.1,0.1,0.1 });
    double polyHTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::VectorXd p = points.col(i);
        Eigen::Vector<double, 4> b0p = b0;
        for (int k = 0; k < n; k++)
        {
            b0p(k) = b0(k) + p(k);
            b0p(k + n) = b0(k + n) - p(k);
        }
        Polyhedron poly(A0, b0p);
        poly.closestPoint(points.col(i+1)); //Compute this to ensure the solver is previously allocated
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            poly.closestPoint(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            polyHTime += ms_double.count();
        }
    }
    std::cout << "PolyhedronH performance" << std::endl;
    std::cout << polyHTime/N << "ms" << std::endl;

    double radius2 = 0.05;
    double obsTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::MatrixXd v = V0.colwise() + points.col(i);
        PolyhedronObstacleCircularRobot obs(v, radius2);
        obs.allocateClosestPointSolver();
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            obs.closestPoint(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            obsTime += ms_double.count();
        }
    }
    std::cout << "PolyhedronObstacleCircularRobot performance" << std::endl;
    std::cout << obsTime/N << "ms" << std::endl;
    std::cout << "All tests passed successful" << std::endl;
}

void isInsidePerformance()
{
    std::cout << "Is Inside Performance" << std::endl;
    /*Performance tests for closest point*/
    int n = 2;
    int N = 100;
    std::chrono::duration<double, std::milli> ms_double;
    Eigen::MatrixXd points = Eigen::MatrixXd::Random(n, N);
    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;
    double pointTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Point p(points.col(i));
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            p.isInside(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            pointTime += ms_double.count();
        }
    }
    std::cout << "Point performance" << std::endl;
    std::cout << pointTime/N << "ms" << std::endl;

    double radius = 0.1;
    double sphereTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Sphere c(points.col(i), radius);
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            c.isInside(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            sphereTime += ms_double.count();
        }
    }
    std::cout << "Sphere performance" << std::endl;
    std::cout << sphereTime/N << "ms" << std::endl;

    Eigen::Matrix<double, 2, 2> C({ {0.1,0.035},{0.035,0.2} });
    Eigen::Vector<double, 2> d({ 0.,0. });
    double ellipsoidTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Ellipsoid e(C, points.col(i));
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            e.isInside(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            ellipsoidTime += ms_double.count();
        }
    }
    std::cout << "Ellipsoid performance" << std::endl;
    std::cout << ellipsoidTime/N << "ms" << std::endl;

    Eigen::Matrix<double, 4, 2> A0({ {1.,0},{0.,1.},{-1.,0},{0.,-1.} });
    Eigen::Vector<double, 4> b0({ 0.1,0.1,0.1,0.1 });
    double polyHTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::VectorXd p = points.col(i);
        Eigen::Vector<double, 4> b0p = b0;
        for (int k = 0; k < n; k++)
        {
            b0p(k) = b0(k) + p(k);
            b0p(k + n) = b0(k + n) - p(k);
        }
        Polyhedron poly(A0, b0p);
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            poly.isInside(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            polyHTime += ms_double.count();
        }
    }
    std::cout << "PolyhedronH performance" << std::endl;
    std::cout << polyHTime/N << "ms" << std::endl;

    Eigen::Matrix<double, 2, 4> V0({ { -0.1, 0.1, 0.1, -0.1}, {-0.1,-0.1,0.1,0.1} });
    double radius2 = 0.05;
    double obsTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::MatrixXd v = V0.colwise() + points.col(i);
        PolyhedronObstacleCircularRobot obs(v, radius2);
        for (int j = i + 1; j < N; j++)
        {
            t1 = std::chrono::high_resolution_clock::now();
            obs.isInside(points.col(j));
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            obsTime += ms_double.count();
        }
    }
    std::cout << "PolyhedronObstacleCircularRobot performance" << std::endl;
    std::cout << obsTime/N << "ms" << std::endl;
    std::cout << "All tests passed successful" << std::endl;
}


void inscribedEllipsoidPerformance()
{
    std::cout << "Inscribed Ellipsoid Performance" << std::endl;
    /*Performance tests for closest point*/
    int n = 2;
    int N = 100;
    std::chrono::duration<double, std::milli> ms_double;
    Eigen::MatrixXd points = Eigen::MatrixXd::Random(n, N);
    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;

    Eigen::Matrix<double, 4, 2> A0({ {1.,0},{0.,1.},{-1.,0},{0.,-1.} });
    Eigen::Vector<double, 4> b0({ 0.1,0.1,0.1,0.1 });
    double polyTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::VectorXd p = points.col(i);
        Eigen::Vector<double, 4> b0p = b0;
        for (int k = 0; k < n; k++)
        {
            b0p(k) = b0(k) + p(k);
            b0p(k + n) = b0(k + n) - p(k);
        }
        Polyhedron poly(A0, b0p);
        t1 = std::chrono::high_resolution_clock::now();
        poly.inscribedEllipsoid();
        t2 = std::chrono::high_resolution_clock::now();
        ms_double = t2 - t1;
        polyTime += ms_double.count();
    }
    std::cout << "Polyhedron performance" << std::endl;
    std::cout << polyTime/N << "ms" << std::endl;
    std::cout << "All tests passed successful" << std::endl;
}

void closestPointExpandingEllipsoidPerformance()
{
    /*Performance tests for closest point expanding ellipsoid*/
    std::cout << "Closest Point Expanding Ellipsoid Performance" << std::endl;
    int n = 2;
    int N = 100;
    std::chrono::duration<double, std::milli> ms_double;
    Eigen::MatrixXd points = Eigen::MatrixXd::Random(n, N);
    Eigen::Matrix<double, 2, 2> C_expanding({ {0.5,0.2},{0.2,0.5}}); //Shape of the expanding ellipsoid
    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;
    double pointTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Point p(points.col(i));
        for (int j = i + 1; j < N; j++)
        {
            Ellipsoid expEllipsoid(C_expanding, points.col(j));
            Eigen::VectorXd p_out(n);
            t1 = std::chrono::high_resolution_clock::now();
            p.closestPointExpandingEllipsoid(expEllipsoid, p_out);
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            pointTime += ms_double.count();
        }
    }
    std::cout << "Point performance" << std::endl;
    std::cout << pointTime / N << "ms" << std::endl;

    double radius = 0.1;
    double sphereTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Sphere c(points.col(i), radius);
        c.allocateClosestPointEllipsoidSolver();
        for (int j = i + 1; j < N; j++)
        {
            Ellipsoid expEllipsoid(C_expanding, points.col(j));
            Eigen::VectorXd p_out(n);
            t1 = std::chrono::high_resolution_clock::now();
            c.closestPointExpandingEllipsoid(expEllipsoid,p_out);
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            sphereTime += ms_double.count();
        }
    }
    std::cout << "Sphere performance" << std::endl;
    std::cout << sphereTime / N << "ms" << std::endl;

    Eigen::Matrix<double, 2, 2> C({ {0.1,0.035},{0.035,0.2} });
    Eigen::Vector<double, 2> d({ 0.,0. });
    double ellipsoidTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Ellipsoid e(C, points.col(i));
        e.allocateClosestPointEllipsoidSolver();
        for (int j = i + 1; j < N; j++)
        {
            Ellipsoid expEllipsoid(C_expanding, points.col(j));
            Eigen::VectorXd p_out(n);
            t1 = std::chrono::high_resolution_clock::now();
            e.closestPointExpandingEllipsoid(expEllipsoid, p_out);
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            ellipsoidTime += ms_double.count();
        }
    }

    std::cout << "Ellipsoid performance" << std::endl;
    std::cout << ellipsoidTime / N << "ms" << std::endl;

    Eigen::Matrix<double, 2, 4> V0({{-0.1, 0.1, 0.1, -0.1}, {-0.1,-0.1,0.1,0.1}});
    double polyVTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::MatrixXd v = V0.colwise() + points.col(i);
        PolyhedronV poly(v);
        poly.allocateClosestPointEllipsoidSolver();
        for (int j = i + 1; j < N; j++)
        {
            Ellipsoid expEllipsoid(C_expanding, points.col(j));
            Eigen::VectorXd p_out(n);
            t1 = std::chrono::high_resolution_clock::now();
            poly.closestPointExpandingEllipsoid(expEllipsoid, p_out);
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            polyVTime += ms_double.count();
        }
    }
    std::cout << "PolyhedronV performance" << std::endl;
    std::cout << polyVTime / N << "ms" << std::endl;

    Eigen::Matrix<double, 4, 2> A0({ {1.,0},{0.,1.},{-1.,0},{0.,-1.} });
    Eigen::Vector<double, 4> b0({ 0.1,0.1,0.1,0.1 });
    double polyHTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::VectorXd p = points.col(i);
        Eigen::Vector<double, 4> b0p = b0;
        for (int k = 0; k < n; k++)
        {
            b0p(k) = b0(k) + p(k);
            b0p(k + n) = b0(k + n) - p(k);
        }
        Polyhedron poly(A0, b0p);
        Ellipsoid expEllipsoid_tmp(C_expanding, points.col(i+1));
        Eigen::VectorXd p_out(n);
        poly.closestPointExpandingEllipsoid(expEllipsoid_tmp, p_out); //Compute this once to ensure the solver is allocated
        for (int j = i + 1; j < N; j++)
        {
            Ellipsoid expEllipsoid(C_expanding, points.col(j));            
            t1 = std::chrono::high_resolution_clock::now();
            poly.closestPointExpandingEllipsoid(expEllipsoid, p_out);
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            polyHTime += ms_double.count();
        }
    }
    std::cout << "PolyhedronH performance" << std::endl;
    std::cout << polyHTime / N << "ms" << std::endl;

    double radius2 = 0.05;
    double obsTime = 0.0;
    for (int i = 0; i < N; i++)
    {
        Eigen::MatrixXd v = V0.colwise() + points.col(i);
        PolyhedronObstacleCircularRobot obs(v, radius2);
        obs.allocateClosestPointEllipsoidSolver();
        for (int j = i + 1; j < N; j++)
        {
            Ellipsoid expEllipsoid(C_expanding, points.col(j));
            Eigen::VectorXd p_out(n);
            t1 = std::chrono::high_resolution_clock::now();
            obs.closestPointExpandingEllipsoid(expEllipsoid, p_out);
            t2 = std::chrono::high_resolution_clock::now();
            ms_double = t2 - t1;
            obsTime += ms_double.count();
        }
    }
    std::cout << "PolyhedronObstacleCircularRobot performance" << std::endl;
    std::cout << obsTime / N << "ms" << std::endl;

    std::cout << "All tests passed successful" << std::endl;
}

int main(int argc, char** argv)
{
    Range* range = Range::getInstance();
    range->setRange(2, -10., 10.);
    //Performance tests for closest point and isInside
    closestPointPerformance();
    isInsidePerformance();
    inscribedEllipsoidPerformance();
    closestPointExpandingEllipsoidPerformance();
    return 0;
}