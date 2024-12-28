#include "CircleNode.h"

CircleNode::CircleNode(const Eigen::VectorXd& center, const double& radius) : circle(center,radius)
{
}

void* CircleNode::getNodeData()
{
	return this;
}
