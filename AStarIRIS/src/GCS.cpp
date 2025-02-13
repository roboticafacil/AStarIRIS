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

void GCS::print(std::ostream& out, const int &num)
{
	int i = 1;
	out << "AGCS" << num << "=cell(0);" << std::endl;
	out << "bGCS" << num << "=cell(0);" << std::endl;
	out << "colGCS" << num << "=cell(0);" << std::endl;
	out << "centroidGCS" << num << "=cell(0);" << std::endl;
	out << "nameGCS" << num << "=cell(0);" << std::endl;
	for (NodeMap::iterator it = this->nodes.begin(); it != this->nodes.end(); it++)
	{
		Node* node = it->second;
		if (dynamic_cast<PolyhedronNode*>(node) != NULL)
		{
			PolyhedronNode* polyNode = (PolyhedronNode*)node;
			polyNode->polyhedron.print(out, "A", "b",true);
			out << "AGCS" << num << "{" << i << "} = -A; " << std::endl;
			out << "bGCS" << num << "{" << i << "}=-b;" << std::endl;
			out << "colGCS" << num << "{" << i << "}=[0.5 0.5 0.5];" << std::endl;
			out << "centroidGCS" << num << "{" << i << "}=mean(con2vert(A,b),2);" << std::endl;
			out << "nameGCS" << num << "{" << i << "}='"<< i <<"';" << std::endl;
			i++;
		}
	}
}

void GCS::print(std::ostream& out, const int& num, const int &agents)
{
	int i = 1;
	out << "AGCS" << num << "=cell(0," << agents << "); " << std::endl;
	out << "bGCS" << num << "=cell(0," << agents << ");" << std::endl;
	out << "colGCS" << num << "=cell(0);" << std::endl;
	out << "centroidGCS" << num << "=cell(0," << agents << ");" << std::endl;
	out << "nameGCS" << num << "=cell(0);" << std::endl;
	for (NodeMap::iterator it = this->nodes.begin(); it != this->nodes.end(); it++)
	{
		Node* node = it->second;
		if (dynamic_cast<PolyhedronNode*>(node) != NULL)
		{
			PolyhedronNode* polyNode = (PolyhedronNode*)node;
			polyNode->polyhedron.print(out, "A", "b", true);
			int n = polyNode->polyhedron.n / agents;
			for (int agent = 1; agent <= agents; agent++)
			{
				out << "p = PolyHedron.projectCoordinatePlaneStatic(A, b," << ((agent - 1) * n + 1) << ":" << agent * n << ");" << std::endl;
				out << "AGCS" << num << "{" << i << "," << agent << "} = -p.A; " << std::endl;
				out << "bGCS" << num << "{" << i << "," << agent << "}=-p.b;" << std::endl;
				out << "centroidGCS" << num << "{" << i << "," << agent << "}=mean(con2vert(p.A,p.b),2);" << std::endl;
			}
			out << "colGCS" << num << "{" << i << "}=[0.5 0.5 0.5];" << std::endl;
			out << "nameGCS" << num << "{" << i << "}='" << i << "';" << std::endl;
			i++;
		}
	}
}

void GCS::printGraph(std::ostream& out, const int& num)
{
	int i = 1;
	out << "g" << num << "=graph();" << std::endl;
	std::map<int, int> nodeMap;
	out << "names" << num << "={";
	for (std::vector<int>::iterator it = this->nodeKeys.begin(); it != this->nodeKeys.end(); it++)
	{
		nodeMap[*it] = i++;
		out << "'" << (*it) + 1 << "' ";
	}
	out << "}';" << std::endl;
	out << "g" << num << "=addnode(g" << num << ",names" << num << ");" << std::endl;
	out << "edges" << num << "=[";
	for (std::vector<Edge>::iterator it = this->edges.begin(); it != this->edges.end(); it++)
	{
		if (it < (this->edges.end() - 2))
			out << nodeMap[it->first] << " " << nodeMap[it->second] << ";" << std::endl;
		else
			out << nodeMap[it->first] << " " << nodeMap[it->second] << std::endl;
	}
	out << "];" << std::endl;
	out << "g" << num << "=addedge(g" << num << ",table(edges" << num << ",'VariableNames',{'EndNodes'}));" << std::endl;
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
