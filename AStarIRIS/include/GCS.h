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
    void print(std::ostream& out, const int &num=1);
    //GCS& operator=(const GCS& other);
    bool contains(const Eigen::VectorXd& q, const double & tol=0.);
    int findConvexSet(const Eigen::VectorXd& q);
    std::vector<int> findConvexSets(const Eigen::VectorXd& q);
    void print();
    void printEdges();
};
#endif