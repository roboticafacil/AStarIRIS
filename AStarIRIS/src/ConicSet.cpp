//#include "fusion.h"
#include "ConicSet.h"

//using namespace mosek::fusion;
//using namespace monty;


ConicSet::ConicSet(const int& n) : n(n)
{
}

ConicSet::ConicSet(const ConicSet& conicSet) : n(conicSet.n), centroidComputed(conicSet.centroidComputed), centroid(conicSet.centroid)
{
}
ConicSet::~ConicSet()
{
}