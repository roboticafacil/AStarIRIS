#include "Point.h"
#include "Node.h"
#include "PointNode.h"

PointNode::PointNode(const Eigen::VectorXd& p) : point(p)
{
}

void* PointNode::getNodeData()
{
	return this;
}