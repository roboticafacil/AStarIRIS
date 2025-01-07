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

void GCS::print(std::ostream& out)
{
	int i = 1;
	for (NodeMap::iterator it = this->nodes.begin(); it != this->nodes.end(); it++)
	{
		Node* node = it->second;
		if (dynamic_cast<PolyhedronNode*>(node) != NULL)
		{
			PolyhedronNode* polyNode = (PolyhedronNode*)node;
			polyNode->polyhedron.print(out, "A", "b",true);
			out << "AGCS{" << i << "}=-A;" << std::endl;
			out << "bGCS{" << i << "}=-b;" << std::endl;
			out << "colGCS{" << i << "}=[0.5 0.5 0.5];" << std::endl;
			out << "centroidGCS{" << i << "}=mean(con2vert(A,b),2);" << std::endl;
			out << "nameGCS{" << i << "}='"<< i <<"';" << std::endl;
			i++;
		}
	}
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

std::vector<int> GCS::findConvexSets(const Eigen::VectorXd& q)
{
	int nodeKey = -1;
	std::vector<int> nodeKeysFound;
	for (std::vector<int>::iterator it = this->nodeKeys.begin(); it != this->nodeKeys.end(); it++)
	{
		Node* node = this->nodes[*it];
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		if (polyNode->polyhedron.isInside(q))
		{
			nodeKeysFound.push_back(*it);
		}
	}
	return nodeKeysFound;
}
