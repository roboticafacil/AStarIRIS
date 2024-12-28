#include "EllipsoidNode.h"

EllipsoidNode::EllipsoidNode(const Eigen::MatrixXd& C, const Eigen::VectorXd& d) : ellipsoid(C,d)
{
}

void* EllipsoidNode::getNodeData()
{
	return this;
}
