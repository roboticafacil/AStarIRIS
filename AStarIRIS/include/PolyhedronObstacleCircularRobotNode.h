#ifndef POLYHEDRON_OBSTACLE_CIRCULAR_ROBOT_NODE
#define POLYHEDRON_OBSTACLE_CIRCULAR_ROBOT_NODE
#include "Node.h"
#include<Eigen/Dense>
#include "PolyhedronObstacleCircularRobot.h"

class PolyhedronObstacleCircularRobotNode: public Node
{
public:
    PolyhedronObstacleCircularRobot polyhedronObstacle;
    PolyhedronObstacleCircularRobotNode(const Eigen::MatrixXd& v, const double& r);
    void* getNodeData();
};
#endif