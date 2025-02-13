#include "NavGraph.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"

//Default constructor
NavGraph::NavGraph() : Graph()
{
}
//Copy constructor
NavGraph::NavGraph(NavGraph* graph) : Graph(graph)
{
}

void NavGraph::printConvexSets(std::ostream& out, const int& num)
{
	int i = 1;
	//out << "ANavGraph" << num << "=cell(0);" << std::endl;
	//out << "bNavGraph" << num << "=cell(0);" << std::endl;
	//out << "colNavGraph" << num << "=cell(0);" << std::endl;
	out << "centroidNavGraph" << num << "=cell(0);" << std::endl;
	out << "nameNavGraph" << num << "=cell(0);" << std::endl;
	std::map<int, int> nodeMap;
	for (std::vector<int>::iterator it = this->nodeKeys.begin(); it != this->nodeKeys.end(); it++)
	{
		nodeMap[*it] = (*it) + 1;
	}
	i = 1;
	for (NodeMap::iterator it = this->nodes.begin(); it != this->nodes.end(); it++)
	{
		Node* node = it->second;
		if (dynamic_cast<PolyhedronNode*>(node) != NULL)
		{
			PolyhedronNode* polyNode = (PolyhedronNode*)node;
			polyNode->polyhedron.print(out, "A", "b", true);
			//if (i == 25)
			//{
			//	std::cout << "Check this out!" << std::endl;
			//	std::cout << polyNode->polyhedron << std::endl;
			//}
			//out << "ANavGraph" << num << "{" << i << "} = -A; " << std::endl;
			//out << "bNavGraph" << num << "{" << i << "}=-b;" << std::endl;
			//out << "colNavGraph" << num << "{" << i << "}=[0.1 0.9 0.1];" << std::endl;
			out << "centroidNavGraph" << num << "{" << i << "}=mean(con2vert(A,b),2);" << std::endl;
			//out << "nameNavGraph" << num << "{" << i << "}='" << nodeMap[it->first] << "';" << std::endl;
			out << "nameNavGraph" << num << "{" << i << "}='" << (it->first+1) << "';" << std::endl;
			
		}
		if (dynamic_cast<PointNode*>(node) != NULL)
		{
			PointNode* pointNode = (PointNode*)node;
			//out << "ANavGraph" << num << "{" << i << "} = []; " << std::endl;
			//out << "bNavGraph" << num << "{" << i << "}=[];" << std::endl;
			//out << "colNavGraph" << num << "{" << i << "}=[0.75 0.25 0.25];" << std::endl;
			pointNode->point.print(out,"p");
			out << "centroidNavGraph" << num << "{" << i << "}=p;" << std::endl;
			//out << "nameNavGraph" << num << "{" << i << "}='" << nodeMap[it->first] << "';" << std::endl;
			out << "nameNavGraph" << num << "{" << i << "}='" << (it->first + 1) << "';" << std::endl;
		}
		i++;
	}
}

void NavGraph::printConvexSet(std::ostream& out, const int& num, const int &nodeKey)
{
	Node* node = this->getNode(nodeKey);
	if (dynamic_cast<PolyhedronNode*>(node) != NULL)
	{
		PolyhedronNode* polyNode = (PolyhedronNode*)node;
		polyNode->polyhedron.print(out, "A", "b", true);
		out << "ANavGraph" << num << " = -A; " << std::endl;
		out << "bNavGraph" << num << "=-b;" << std::endl;
		out << "colNavGraph" << num << "=[0.1 0.9 0.1];" << std::endl;
	}
	else if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
	{
		PolyhedronNode* polyNode = (PolyhedronNode*)node;
		polyNode->polyhedron.print(out, "A", "b", true);
		out << "ANavGraph" << num << " = -A; " << std::endl;
		out << "bNavGraph" << num << "=-b;" << std::endl;
		out << "colNavGraph" << num << "=[0.1 0.9 0.1];" << std::endl;
	}
	else if (dynamic_cast<PointNode*>(node) != NULL)
	{
		PointNode* pointNode = (PointNode*)node;
		out << "ANavGraph" << num << " = []; " << std::endl;
		out << "bNavGraph" << num << " =[];" << std::endl;
		out << "colNavGraph" << num << " =[0.75 0.25 0.25];" << std::endl;
	}
}

void NavGraph::printConvexSet(std::ostream& out, const int& num, const int& nodeKey, const int &agents)
{
	Node* node = this->getNode(nodeKey);
	if (dynamic_cast<PolyhedronNode*>(node) != NULL)
	{
		PolyhedronNode* polyNode = (PolyhedronNode*)node;
		int n = polyNode->polyhedron.n/agents;
		polyNode->polyhedron.print(out, "A", "b", true);
		for (int agent = 1; agent <= agents; agent++)
		{
			out << "p = PolyHedron.projectCoordinatePlaneStatic(A, b," << ((agent - 1) * n + 1) << ":" << agent * n << ");" << std::endl;
			out << "ANavGraph" << num << "{" << agent << "} = -p.A; " << std::endl;
			out << "bNavGraph" << num << "{" << agent << "} =-p.b;" << std::endl;
		}
		out << "colNavGraph" << num << "=[0.1 0.9 0.1];" << std::endl;
	}
	else if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
	{
		PolyhedronNode* polyNode = (PolyhedronNode*)node;
		polyNode->polyhedron.print(out, "A", "b", true);
		for (int agent = 1; agent <= agents; agent++)
		{
			out << "ANavGraph" << num << "{" << agent << "} = -A; " << std::endl;
			out << "bNavGraph" << num << "{" << agent << "}=-b;" << std::endl;
		}
		out << "colNavGraph" << num << "=[0.1 0.9 0.1];" << std::endl;
	}
	else if (dynamic_cast<PointNode*>(node) != NULL)
	{
		PointNode* pointNode = (PointNode*)node;
		for (int agent = 1; agent <= agents; agent++)
		{
			out << "ANavGraph" << num << "{" << agent << "} = []; " << std::endl;
			out << "bNavGraph" << num << "{" << agent << "} =[];" << std::endl;
		}
		out << "colNavGraph" << num << " =[0.75 0.25 0.25];" << std::endl;
	}
}

void NavGraph::printGraph(std::ostream& out, const int& num,std::vector<double> &weights)
{
	int i = 1;
	out << "g" << num << "=digraph();" << std::endl;
	std::map<int,int> nodeMap;
	out << "names" << num << "={";
	for (std::vector<int>::iterator it = this->nodeKeys.begin(); it != this->nodeKeys.end(); it++)
	{
		nodeMap[*it] = i++;
		out << "'" << (*it)+1 << "' ";
	}
	out << "}';" << std::endl;
	out << "g" << num << "=addnode(g" << num << ",names" << num << ");" << std::endl;
	out << "edges" << num << "=[";
	for (std::vector<Edge>::iterator it = this->edges.begin(); it != this->edges.end(); it++)
	{
		if (it<(this->edges.end() - 2))
			out << nodeMap[it->first] << " " << nodeMap[it->second] << ";" << std::endl;
		else
			out << nodeMap[it->first] << " " << nodeMap[it->second] << std::endl;
	}
	out << "];" << std::endl;
	if (weights.size() == this->edges.size())
	{
		out << "weights" << num << "=[";
		for (std::vector<double>::iterator it = weights.begin(); it != weights.end(); it++)
		{
			if (it < (weights.end() - 2))
				out << (*it) << ";" << std::endl;
			else
				out << (*it) << std::endl;
		}
		out << "];" << std::endl;
	}
	if (weights.size() == this->edges.size())
		out << "g" << num << "=addedge(g" << num << ",table(edges" << num << ",weights"<<num<< ", 'VariableNames', {'EndNodes','Weights'})); " << std::endl;
	else
		out << "g" << num << "=addedge(g" << num << ",table(edges" << num << ",'VariableNames',{'EndNodes'}));" << std::endl;
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

void NavGraph::print()
{
	Graph::print();
}

bool NavGraph::contains(const Eigen::VectorXd& q, const double& tol)
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
