#include "EGCS.h"
#include "PolyhedronNode.h"

//Default constructor
EGCS::EGCS() : Graph()
{
}
//Copy constructor
EGCS::EGCS(EGCS* graph) : Graph(graph)
{
}

/*void EGCS::buildEdgeGraph(GCS& gcs)
{
	for (int i = 0; i < gcs.numEdges; i++)
	{
		Edge edge = gcs.getEdge(i);
		PolyhedronNode* nodeA = (PolyhedronNode*)gcs.getNode(edge.first);
		PolyhedronNode* nodeB = (PolyhedronNode*)gcs.getNode(edge.second);
		//Vertically concatenate the polyhedra constraints
		Eigen::MatrixXd bigA(nodeA->polyhedron.A.rows() + nodeB->polyhedron.A.rows(), nodeA->polyhedron.A.cols());
		Eigen::VectorXd bigB(nodeA->polyhedron.A.rows() + nodeB->polyhedron.A.rows());
		bigA << nodeA->polyhedron.A, nodeB->polyhedron.A;
		bigB << nodeA->polyhedron.b, nodeB->polyhedron.b;
		PolyhedronNode* node = new PolyhedronNode(bigA, bigB);
		this->addNode(node);
		//nodeKeyPairs.push_back(std::make_pair(edge.first, edge.second));
	}
	//std::cout << "Debugging... " << std::endl;
	std::map<std::pair<int, int>, int> edgeMap;
	for (int i = 0; i < gcs.numEdges; i++)
	{
		Edge edge1 = gcs.getEdge(i);
		for (int j = 0; j < gcs.numEdges; j++)
		{
			if (i != j)
			{
				Edge edge2 = gcs.getEdge(j);
				if (edge2.first == edge1.first)
				{
					//std::stringstream s;
					//s << nodeKeys[i] << "->" << nodeKeys[j];
					edgeMap.emplace(std::make_pair(nodeKeys[i], nodeKeys[j]), nodeKeys[i]);
					//edgeMap.emplace(std::make_pair(nodeKeys[i], edge1.first), nodeKeys[j]);
					//this->addEdge(nodeKeys[i], nodeKeys[j]);
					continue;
				}
				if (edge2.second == edge1.first)
				{
					//std::stringstream s;
					//s << nodeKeys[i] << "->" << nodeKeys[j];
					edgeMap.emplace(std::make_pair(nodeKeys[i], nodeKeys[j]), nodeKeys[i]);
					//edgeMap.emplace(std::make_pair(nodeKeys[i], edge1.first), nodeKeys[j]);
					//this->addEdge(nodeKeys[i], nodeKeys[j]);
					continue;
				}
			}
		}
		for (int j = 0; j < gcs.numEdges; j++)
		{
			if (i != j)
			{
				Edge edge2 = gcs.getEdge(j);
				if (edge2.first == edge1.second)
				{
					//std::stringstream s;
					//s << nodeKeys[i] << "->" << nodeKeys[j];
					edgeMap.emplace(std::make_pair(nodeKeys[i], nodeKeys[j]), nodeKeys[i]);
					//edgeMap.emplace(std::make_pair(nodeKeys[i], edge1.second), nodeKeys[j]);
					//this->addEdge(nodeKeys[i], nodeKeys[j]);
					continue;
				}
				if (edge2.second == edge1.second)
				{
					//std::stringstream s;
					//s << nodeKeys[i] << "->" << nodeKeys[j];
					edgeMap.emplace(std::make_pair(nodeKeys[i], nodeKeys[j]), nodeKeys[i]);
					//edgeMap.emplace(std::make_pair(nodeKeys[i], edge1.second), nodeKeys[j]);
					//this->addEdge(nodeKeys[i], nodeKeys[j]);
					continue;
				}
			}
		}
	}

	for (std::map<std::pair<int, int>, int>::iterator it = edgeMap.begin(); it != edgeMap.end(); it++)
	{
		this->addEdge(it->first.first, it->first.second);
	}
}*/

/*EGCS& EGCS::operator=(const EGCS& other)
{
	*this=other;
	return *this;
}*/

bool EGCS::contains(const Eigen::VectorXd& q, const double& tol)
{
	for (std::vector<int>::iterator it = this->nodeKeys.begin(); it != this->nodeKeys.end(); it++)
	{
		Node* node = this->nodes[*it];
		PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
		if (polyNode->polyhedron.isInside(q, tol))
			return true;
	}
	return false;
}
