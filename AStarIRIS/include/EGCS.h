#ifndef EGCS_H
#define EGCS_H
#include "Graph.h"
#include "GCS.h"

class EGCS : public Graph
{
public:
    //Default constructor
    EGCS();
    //Copy constructor
    EGCS(EGCS* graph);
    //EGCS& operator=(const EGCS& other);
    //void buildEdgeGraph(GCS& gcs);
    bool contains(const Eigen::VectorXd& q, const double& tol = 0.);
};
#endif