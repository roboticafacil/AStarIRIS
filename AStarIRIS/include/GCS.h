#ifndef GCS_H
#define GCS_H
#include "Graph.h"

class GCS : public Graph
{
public:
    //Default constructor
    GCS();
    //Copy constructor
    GCS(GCS* graph);
    //GCS constructor
    void buildEdgeGraph(Graph* graph);
    bool contains(const Eigen::VectorXd& q);
};
#endif