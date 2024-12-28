#ifndef ELLIPSOID_NODE
#define ELLIPSOID_NODE
#include<Eigen/Dense>
#include "Ellipsoid.h"
#include "Node.h"
class EllipsoidNode: public Node
{
public:
	Ellipsoid ellipsoid;
	EllipsoidNode(const Eigen::MatrixXd& C, const Eigen::VectorXd& d);
	void* getNodeData();
};
#endif