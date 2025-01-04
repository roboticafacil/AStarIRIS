#ifndef POLYHEDRON_NODE
#define POLYHEDRON_NODE
#include<Eigen/Dense>
#include "Polyhedron.h"
#include "Node.h"

class PolyhedronNode : public Node
{
public:
    Polyhedron polyhedron;
    Polyhedron shrinkedPolyhedron;
    bool useShrinked;
    //std::vector<bool> facetExpanded;
    PolyhedronNode(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
    PolyhedronNode(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& bShrinked);
    void* getNodeData();
};
#endif