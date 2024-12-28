#include<Eigen/Dense>
#include "Polyhedron.h"
#include "Node.h"
#include "PolyhedronNode.h"

PolyhedronNode::PolyhedronNode(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& bShrinked) : polyhedron(A, b), shrinkedPolyhedron(A,bShrinked), useShrinked(true)
{
}

PolyhedronNode::PolyhedronNode(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) : polyhedron(A, b), shrinkedPolyhedron(A.cols()), useShrinked(false)
{
}

void* PolyhedronNode::getNodeData() {
	return this;
}