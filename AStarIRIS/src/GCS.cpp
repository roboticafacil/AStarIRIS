#include "GCS.h"
#include "PolyhedronNode.h"

//Default constructor
GCS::GCS() : Graph()
{
}
//Copy constructor
GCS::GCS(GCS* graph): Graph(graph)
{
}

/*GCS& GCS::operator=(const GCS& other)
{
	*this = other;
	return *this;
}*/

void GCS::print()
{
	Graph::printNodes();
	this->printEdges();
}

void GCS::printEdges()
{
	std::cout << "Graph edges: " << this->numEdges << std::endl;
	for (std::vector<Edge>::iterator it = this->edges.begin(); it != this->edges.end(); ++it)
	{
		//std::cout << "edge: " << std::endl;
		std::cout << "edge: " << it->first << "<->" << it->second << std::endl;
	}
}

bool GCS::contains(const Eigen::VectorXd& q, const double &tol)
{
	for (std::vector<int>::iterator it = this->nodeKeys.begin(); it!=this->nodeKeys.end(); it++)
	{
		Node* node = this->nodes[*it];
		PolyhedronNode* polyNode=(PolyhedronNode*)node->getNodeData();
		if (polyNode->polyhedron.isInside(q,tol))
			return true;
	}
	return false;
}

int GCS::findConvexSet(const Eigen::VectorXd& q)
{
	int nodeKey = -1;
	for (std::vector<int>::iterator it = this->nodeKeys.begin(); it != this->nodeKeys.end(); it++)
	{
		Node* node = this->nodes[*it];
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		if (polyNode->polyhedron.isInside(q))
		{
			nodeKey = *it;
			break;
		}
	}
	return nodeKey;
}
