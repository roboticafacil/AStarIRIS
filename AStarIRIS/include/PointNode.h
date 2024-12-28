#ifndef POINT_NODE
#define POINT_NODE
#include<Eigen/Dense>
#include "Point.h"
#include "Node.h"

class PointNode : public Node
{
public:
    Point point;
    PointNode(const Eigen::VectorXd& p);
    void* getNodeData();
};
#endif