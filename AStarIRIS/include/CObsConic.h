#ifndef COBS_CONIC
#define COBS_CONIC
#include "CObs.h"
#include <list>
#include "ConicSet.h"
#include "Point.h"
#include "Circle.h"
#include "Polyhedron.h"
#include "PolyhedronObstacleCircularRobot.h"

class CObsConic: public CObs
{
public:
	CObsConic();
	CObsConic(CObsConic* cObs);
	void addObject(ConicSet* object);
	bool isFree(const Eigen::VectorXd& q);
	/*void addObject(Point* point);
	void addObject(Circle* circle);
	void addObject(Ellipsoid* ellipsoid);
	void addObject(Polyhedron* polyhedron);
	void addObject(PolyhedronObstacleCircularRobot* polyhedronCircularRobot);*/
	int numObjects=0;
	/*std::vector<Point*> points;
	std::vector<Circle*> circles;
	std::vector<Ellipsoid*> ellipsoids;
	std::vector<Polyhedron*> polyhedra;
	std::vector<PolyhedronObstacleCircularRobot*> polyhedraCircularRobot;*/
	std::vector<ConicSet*> conicObjects;	
};
#endif