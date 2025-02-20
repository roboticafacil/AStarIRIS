#ifndef GRAPH
#define GRAPH
#include<Eigen/Dense>
#include <limits>
#include <list>
#include <vector>
#include <map>
#include "Node.h"
#include "Edge.h"

typedef std::map<int, Node*> NodeMap;
typedef std::pair<int, int> NodePair;
typedef std::map<NodePair, Edge*> EdgeMap;

class Graph
{
public:
    int numNodes;
    int numEdges;
protected:
    NodeMap nodes;
    EdgeMap edges;
    std::vector<NodePair> nodePairs;
    std::vector<int> nodeKeys;
public:
    //Default constructor
    Graph();
    //Copy constructor
    Graph(Graph* graph);
    Graph& operator=(Graph& other);
    virtual void print();
    virtual void printNodes();
    virtual void printEdges(); 
    Graph subGraph(std::vector<int>& keys);
    std::vector<int> getNodeKeys();
    Node* getNode(const int& key);
    std::vector<Node*> getNodes();
    std::vector<Node*> getNodes(std::vector<int>& keys);
    int addNode(const Node* node=NULL);
    void addNode(const int& key, const Node* node=NULL);
    bool removeNode(const int& key);

    NodePair addEdge(const int& keyFrom, const int& keyTo, const Edge* edge=NULL);
    void addEdge(const NodePair& nodePair, const Edge* edge=NULL);
    //bool addEdge(const int& keyFrom, const int& keyTo);
    bool removeEdge(const NodePair& nodePair);
    bool removeEdge(const int& keyFrom, const int& keyTo);
    Edge* getEdge(const NodePair& nodePair);
    NodePair getEdgeNodePair(const int& edgeIdx);
    Edge* getEdge(const int& edgeIdx);
    std::vector<NodePair> getEdgeNodePairs();
    
    std::vector<int> findOutEdges(const int& key);
    std::vector<int> findInEdges(const int& key);
    std::vector<Node*> getOutNeighbours(const int& key);
    std::vector<Node*> getInNeighbours(const int& key);
    std::vector<Node*> getNeighbours(const int& key);
    std::vector<int> getOutNeighbourKeys(const int& key);
    std::vector<int> getInNeighbourKeys(const int& key);
    std::vector<int> getNeighbourKeys(const int& key);
    //NodePair getEdge(const int &edgeIdx);
    void setEdges(const std::vector<NodePair>& edges);
};
#endif