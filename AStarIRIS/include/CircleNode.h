#ifndef CIRCLE_NODE
#define CIRCLE_NODE
#include "Node.h"
#include "Circle.h"
class CircleNode :   public Node
{
public:
	Circle circle;
	CircleNode(const Eigen::VectorXd& center, const double &radius);
	void* getNodeData();
};
#endif