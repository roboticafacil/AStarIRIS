#include<Eigen/Dense>
#include "Polyhedron.h"
#include "Node.h"
#include "PolyhedronVNode.h"

PolyhedronVNode::PolyhedronVNode(const Eigen::MatrixXd& v) : polyhedron(v)
{
}

void* PolyhedronVNode::getNodeData() {
	return &this->polyhedron;
}