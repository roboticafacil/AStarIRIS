#ifndef GRAPH
#define GRAPH
#include<Eigen/Dense>
#include <limits>
#include <list>
#include <vector>
#include <map>
#include "Node.h"

typedef std::map<int, Node*> NodeMap;
typedef std::pair<int, int> Edge;

class Graph
{
public:
    int numNodes;
    int numEdges;
protected:
    NodeMap nodes;
    std::vector<Edge> edges;
    std::vector<int> nodeKeys;
public:
    //Default constructor
    Graph();
    //Copy constructor
    Graph(Graph* graph);
    void print();
    void printNodes();
    void printEdges();
    std::vector<int> getNodeKeys();
    Node* getNode(const int& key);
    std::vector<Node*> getNodes();
    std::vector<Node*> getNodes(std::vector<int>& keys);
    Graph subGraph(std::vector<int>& keys);
    int addNode(const Node* node);
    void addNode(const Node* node, const int& key);
    bool removeNode(const int& key);
    bool addEdge(const int& keyFrom, const int& keyTo);
    bool removeEdge(const int& keyFrom, const int& keyTo);
    std::vector<int> findOutEdges(const int& key);
    std::vector<int> findInEdges(const int& key);
    std::vector<Edge> getEdges();
    std::vector<Node*> getOutNeighbours(const int& key);
    std::vector<Node*> getInNeighbours(const int& key);
    std::vector<Node*> getNeighbours(const int& key);
    std::vector<int> getOutNeighbourKeys(const int& key);
    std::vector<int> getInNeighbourKeys(const int& key);
    std::vector<int> getNeighbourKeys(const int& key);
    Edge getEdge(const int &edgeIdx);
};
#endif