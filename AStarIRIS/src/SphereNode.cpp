#include "SphereNode.h"

SphereNode::SphereNode(const Eigen::VectorXd& center, const double& radius) : sphere(center,radius)
{
}

void* SphereNode::getNodeData()
{
	return this;
}
