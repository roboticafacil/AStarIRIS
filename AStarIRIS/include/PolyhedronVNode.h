#ifndef POLYHEDRONV_NODE
#define POLYHEDRONV_NODE
#include<Eigen/Dense>
#include "PolyhedronV.h"
#include "Node.h"

class PolyhedronVNode : public Node
{
public:
    PolyhedronV polyhedron;
    PolyhedronVNode(const Eigen::MatrixXd& v);
    void* getNodeData();
};
#endif