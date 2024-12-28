//#include "fusion.h"
#include "ConicSet.h"

//using namespace mosek::fusion;
//using namespace monty;


ConicSet::ConicSet(const int& n) : n(n)
{
}


/*ConicSet::ConicSet(const int& n, const bool& initializeModel) : n(n), initializeModel(initializeModel)
{
    if (initializeModel)
    {
        M = new Model("ConicSet");
        x = M->variable("x", n, Domain::unbounded());
    }
    else
        M = nullptr;
    //
}*/

//ConicSet::ConicSet(const ConicSet& conicSet): n(conicSet.n), initializeModel(conicSet.initializeModel), centroidComputed(conicSet.centroidComputed), centroid(conicSet.centroid)
ConicSet::ConicSet(const ConicSet& conicSet) : n(conicSet.n), centroidComputed(conicSet.centroidComputed), centroid(conicSet.centroid)
{
}
ConicSet::~ConicSet()
{
    //if (initializeModel)
    //    auto _M = finally([&]() { M->dispose(); });
}