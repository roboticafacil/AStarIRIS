#ifndef SPHERE_NODE
#define SPHERE_NODE
#include "Node.h"
#include "Sphere.h"
class SphereNode :   public Node
{
public:
	Sphere sphere;
	SphereNode(const Eigen::VectorXd& center, const double &radius);
	void* getNodeData();
};
#endif