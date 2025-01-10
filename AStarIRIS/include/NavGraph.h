#ifndef EGCS_H
#define EGCS_H
#include "Graph.h"
#include "GCS.h"

class NavGraph : public Graph
{
public:
    //Default constructor
    NavGraph();
    //Copy constructor
    NavGraph(NavGraph* graph);
    void printConvexSets(std::ostream& out, const int& num = 1);
    void printConvexSet(std::ostream& out, const int& num, const int& nodeKey);
    void printGraph(std::ostream& out, const int& num = 1, std::vector<double>& weights=std::vector<double>());
    //NavGraph& operator=(const NavGraph& other);
    //void buildEdgeGraph(GCS& gcs);
    bool contains(const Eigen::VectorXd& q, const double& tol = 0.);
    void print();
};
#endif