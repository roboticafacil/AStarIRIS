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

void GCS::buildEdgeGraph(Graph* graph)
{
	for (int i = 0; i < graph->numEdges; i++)
	{
		Edge edge=graph->getEdge(i);
		PolyhedronNode* nodeA = (PolyhedronNode*)graph->getNode(edge.first);
		PolyhedronNode* nodeB = (PolyhedronNode*)graph->getNode(edge.second);
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
	std::map<std::pair<int,int>,int> edgeMap;
	for (int i = 0; i < graph->numEdges; i++)
	{
		Edge edge1 = graph->getEdge(i);
		for (int j = 0; j < graph->numEdges; j++)
		{
			if (i!=j)
			{ 
				Edge edge2 = graph->getEdge(j);
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
		for (int j = 0; j < graph->numEdges; j++)
		{
			if (i != j)
			{
				Edge edge2 = graph->getEdge(j);
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
//		std::cout << it->first.first << "->"<< it->first.second << std::endl;
	}
	/*for (int i = 0; i < graph->numNodes; i++)
	{
		std::vector<int> neighbourKeys=graph->getNeighbourKeys(keys[i]);
		for (int j = 0; j < neighbourKeys.size(); j++)
		{
			for (int k = 0; k < graph->numEdges; k++)
			{
				Edge edge1 = graph->getEdge(k);
				if (edge1.first == neighbourKeys[j] && edge1.second == keys[i])
				{
					for (int l = 0; l < neighbourKeys.size(); l++)
					{
						for (int m = 0; m < graph->numEdges; m++)
						{
							Edge edge2 = graph->getEdge(k);
							if (edge2.first == neighbourKeys[l] && edge2.second == keys[i])
							{
								this->addEdge(nodeKeys[k], nodeKeys[m]);
								break;
							}
							if (edge2.second == neighbourKeys[l] && edge2.first == keys[i])
							{
								this->addEdge(nodeKeys[k], nodeKeys[m]);
								break;
							}
						}
					}
					break;
				}
				if (edge1.second == neighbourKeys[j] && edge1.first == keys[i])
				{
					for (int l = 0; l < neighbourKeys.size(); l++)
					{
						for (int m = 0; m < graph->numEdges; m++)
						{
							Edge edge2 = graph->getEdge(k);
							if (edge2.first == neighbourKeys[l] && edge2.second == keys[i])
							{
								this->addEdge(nodeKeys[k], nodeKeys[m]);
								break;
							}
							if (edge2.second == neighbourKeys[l] && edge2.first == keys[i])
							{
								this->addEdge(nodeKeys[k], nodeKeys[m]);
								break;
							}
						}
					}
					break;
				}
			}
		}
	}*/
}

bool GCS::contains(const Eigen::VectorXd& q)
{
	for (std::vector<int>::iterator it = this->nodeKeys.begin(); it!=this->nodeKeys.end(); it++)
	{
		Node* node = this->nodes[*it];
		PolyhedronNode* polyNode=(PolyhedronNode*)node->getNodeData();
		if (polyNode->polyhedron.isInside(q))
			return true;
	}
	return false;
}
