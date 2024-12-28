#include "CObsConic.h"

CObsConic::CObsConic()
{
}

CObsConic::CObsConic(CObsConic* cObs)
{
	this->conicObjects = cObs->conicObjects;
	this->numObjects = cObs->numObjects;
}

void CObsConic::addObject(ConicSet* object)
{
	conicObjects.push_back(object);
	numObjects++;
}

bool CObsConic::isFree(const Eigen::VectorXd& q)
{
	for (std::vector<ConicSet*>::iterator it = conicObjects.begin(); it != conicObjects.end(); it++)
	{
		ConicSet* conicSet = *it;
		if (conicSet->isInside(q))
			return false;
	}
	return true;
}

/*
void CObsConic::addObject(Point* point)
{
	points.push_back(point);
	this->addObject(point);
}

void CObsConic::addObject(Circle* circle)
{
	circles.push_back(circle);
	this->addObject(circle);
}

void CObsConic::addObject(Ellipsoid* ellipsoid)
{
	ellipsoids.push_back(ellipsoid);
	this->addObject(ellipsoid);
}

void CObsConic::addObject(Polyhedron* polyhedron)
{
	polyhedra.push_back(polyhedron);
	this->addObject(polyhedron);
}

void CObsConic::addObject(PolyhedronObstacleCircularRobot* polyhedronCircularRobot)
{
	polyhedraCircularRobot.push_back(polyhedronCircularRobot);
	this->addObject(polyhedronCircularRobot);
}*/
