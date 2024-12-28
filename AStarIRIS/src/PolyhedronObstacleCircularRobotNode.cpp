#include "PolyhedronObstacleCircularRobotNode.h"

PolyhedronObstacleCircularRobotNode::PolyhedronObstacleCircularRobotNode(const Eigen::MatrixXd& v, const double& r): polyhedronObstacle(v,r)
{

}


void* PolyhedronObstacleCircularRobotNode::getNodeData()
{
	return this;
}

