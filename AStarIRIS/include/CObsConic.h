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
	bool isFree(const Eigen::VectorXd& q, const double &tol);
	void print(std::ostream& out);
	int numObjects=0;
	std::vector<ConicSet*> conicObjects;	
};
#endif