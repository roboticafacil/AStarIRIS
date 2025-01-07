#include <iostream>
#include<Eigen/Dense>
#include <limits>
#include <list>
#include <vector>
#include <map>
#include "Graph.h"
#include "Node.h"


    //Default constructor
    Graph::Graph()
    {
        this->nodes = NodeMap();
        this->edges = std::vector<Edge>(0);
        this->numEdges = 0;
        this->numNodes = 0;
    };
    //Copy constructor
    Graph::Graph(Graph* graph) : numEdges(graph->numEdges), numNodes(graph->numNodes), nodes(graph->nodes), edges(graph->edges)
    {
    }

    Graph& Graph::operator=(Graph& other)
    {
        this->numEdges = other.numEdges;
        this->numNodes = other.numNodes;
        this->nodes.clear();
        NodeMap::iterator it = this->nodes.begin();
        for (NodeMap::iterator itOther = other.nodes.begin(); itOther != other.nodes.end(); ++itOther)
        {
            this->nodes[itOther->first] = itOther->second;
        }
        this->nodes = other.nodes;
        this->edges = other.edges;
        this->nodeKeys = other.nodeKeys;
        return *this;
    }

    void Graph::print()
    {
        printNodes();
        printEdges();
    }

    void Graph::printNodes()
    {
        std::cout << "Graph nodes: " << this->numNodes << std::endl;
        for (NodeMap::iterator it = this->nodes.begin(); it != this->nodes.end(); ++it)
        {
            //std::cout << "key: " << std::endl;
            std::cout << "key: " << it->first << std::endl;
        }
    }

    void Graph::printEdges()
    {
        std::cout << "Graph edges: " << this->numEdges << std::endl;
        for (std::vector<Edge>::iterator it = this->edges.begin(); it != this->edges.end(); ++it)
        {
            //std::cout << "edge: " << std::endl;
            std::cout << "edge: " << it->first << "->" << it->second << std::endl;
        }
    }
    std::vector<int> Graph::getNodeKeys()
    {
        return nodeKeys;
    }

    Node* Graph::getNode(const int& key)
    {
        /*NodeMap::iterator it = this->nodes.find(key);
        if (it != this->nodes.end())
        {
            return it->second;
        }
        return NULL;*/
        return this->nodes[key];
    }

    std::vector<Node*> Graph::getNodes()
    {
        std::vector<Node*> _nodes;
        std::vector<int>::iterator it = nodeKeys.begin();
        for (; it != nodeKeys.end(); it++)
        {
            //nodes.push_back(this->getNode(*it));
            _nodes.push_back(this->nodes[*it]);
        }
        return _nodes;
    }

    std::vector<Node*> Graph::getNodes(std::vector<int>& keys)
    {
        std::vector<Node*> nodes;
        std::vector<int>::iterator it = keys.begin();
        for (; it != keys.end(); it++)
        {
            nodes.push_back(this->getNode(*it));
        }
        return nodes;
    }

    Graph Graph::subGraph(std::vector<int>& keys)
    {
        std::vector<Node*> nodes = this->getNodes(keys);
        Graph _graph;
        std::vector<int> _nodeKeys;  //No se usa??
        std::vector<int> outEdges;
        std::vector<int> inEdges;
        for (int i=0;i<keys.size();i++)
        {
            _graph.addNode(nodes[i], keys[i]);
        }
        for (int i = 0; i < keys.size(); i++)
        {
            std::vector<int> oE=this->findOutEdges(keys[i]);
            for (int j = 0; j < oE.size(); j++)
            {
                Edge e=this->getEdge(oE[j]);
                _graph.addEdge(e.first, e.second);
            }
            /*std::vector<int> iE = this->findInEdges(nodeKeys[i]);
            for (int j = 0; j < iE.size(); j++)
            {
                Edge e = this->getEdge(iE[j]);
                _graph.addEdge(e.first, e.second);
            }*/
        }
        return _graph;
    }

    int Graph::addNode(const Node* node)
    {
        //int key = this->nodes.size();
        int key;
        if (this->nodeKeys.size() > 0)
        {
            std::vector<int>::iterator it = this->nodeKeys.end() - 1;
            key = (*it) + 1;
        }
        else
            key = 0;
        this->nodes[key] = (Node*)node;
        this->nodeKeys.push_back(key);
        this->numNodes++;
        return key;
    };

    void Graph::addNode(const Node* node, const int& key)
    {
        this->nodes[key] = (Node*)node;
        this->nodeKeys.push_back(key);
        this->numNodes++;
    };



    bool Graph::removeNode(const int& key)
    {
        NodeMap::iterator it = this->nodes.find(key);
        if (it != this->nodes.end())
        {
            //Remove key from nodeKeys
            std::vector<int>::iterator keyIt = this->nodeKeys.begin();
            for (; keyIt != this->nodeKeys.end(); keyIt++)
            {
                if (*keyIt == key)
                {
                    this->nodeKeys.erase(keyIt);
                    break;
                }
            }
            //First removes the edges that contain that node
            std::vector<Edge>::iterator edgeIt = this->edges.begin();
            while (edgeIt != this->edges.end())
            {
                if (edgeIt->first == key)
                {
                    this->edges.erase(edgeIt);
                    this->numEdges--;
                    continue;
                }
                if (edgeIt->second == key)
                {
                    this->edges.erase(edgeIt);
                    this->numEdges--;
                    continue;
                }
                ++edgeIt;
            }
            //Now deletes that node
            this->nodes.erase(it);
            this->numNodes--;
        }
        return false;
    }

    bool Graph::addEdge(const int& keyFrom, const int& keyTo)
    {
        NodeMap::iterator itFrom = this->nodes.find(keyFrom);
        NodeMap::iterator itTo = this->nodes.find(keyTo);
        if ((itFrom != this->nodes.end()) && (itTo != this->nodes.end()))
        {
            //Avoid adding duplicates
            for (std::vector<std::pair<int, int>>::iterator it = this->edges.begin(); it != this->edges.end(); it++)
            {
                if ((it->first == keyFrom) && (it->second == keyTo))
                    return false;
            }
            Edge edge = std::make_pair(keyFrom, keyTo);
            this->edges.push_back(edge);
            this->numEdges++;
            return true;
        }
        return false;
    };

    bool Graph::removeEdge(const int& keyFrom, const int& keyTo)
    {
        NodeMap::iterator itFrom = this->nodes.find(keyFrom);
        NodeMap::iterator itTo = this->nodes.find(keyTo);
        if ((itFrom != this->nodes.end()) && (itTo != this->nodes.end()))
        {
            for (std::vector<Edge>::iterator edgeIt = this->edges.begin(); edgeIt != this->edges.end(); ++edgeIt)
            {
                if ((edgeIt->first == keyFrom) && (edgeIt->second == keyTo))
                {
                    this->edges.erase(edgeIt);
                    this->numEdges--;
                    return true;
                }
            }
        }
        return false;
    }

    std::vector<int> Graph::findOutEdges(const int& key)
    {
        std::vector<int> _edges;
        for (std::vector<Edge>::iterator edgeIt = this->edges.begin(); edgeIt != this->edges.end(); ++edgeIt)
        {
            if (edgeIt->first == key)
            {
                _edges.push_back(edgeIt- this->edges.begin());
            }
        }
        return _edges;
    }


    std::vector<int> Graph::findInEdges(const int& key)
    {
        std::vector<int> _edges;
        for (std::vector<Edge>::iterator edgeIt = this->edges.begin(); edgeIt != this->edges.end(); ++edgeIt)
        {
            if (edgeIt->second == key)
            {
                _edges.push_back(edgeIt- this->edges.begin());
            }
        }
        return _edges;
    }

    std::vector<Edge> Graph::getEdges()
    {
        return this->edges;
    }

    Edge Graph::getEdge(const int& edgeIdx)
    {
        if (edgeIdx < this->edges.size() && edgeIdx >= 0)
        {
            return this->edges[edgeIdx];
        }
        else
            return Edge(-1,1);
    }

    void Graph::setEdges(const std::vector<Edge>& edges)
    {
        this->edges = edges;
        this->numEdges = this->edges.size();
    }

    std::vector<Node*> Graph::getOutNeighbours(const int& key)
    {
        std::vector<Node*> _nodes;
        std::vector<int> _edges = this->findOutEdges(key);
        for (std::vector<int>::iterator edgeIt = _edges.begin(); edgeIt != _edges.end(); ++edgeIt)
        {
            _nodes.push_back(this->getNode(this->edges[*edgeIt].second));
        }
        return _nodes;
    }

    std::vector<Node*> Graph::getInNeighbours(const int& key)
    {
        std::vector<Node*> _nodes;
        std::vector<int> _edges = this->findInEdges(key);
        for (std::vector<int>::iterator edgeIt = _edges.begin(); edgeIt != _edges.end(); ++edgeIt)
        {
            _nodes.push_back(this->getNode(this->edges[*edgeIt].first));
        }
        return _nodes;
    }

    std::vector<Node*> Graph::getNeighbours(const int& key)
    {
        std::vector<Node*> _nodes;
        std::vector<int> _edgesIn = this->findInEdges(key);
        std::vector<int> _edgesOut = this->findOutEdges(key);
        for (std::vector<int>::iterator edgeIt = _edgesIn.begin(); edgeIt != _edgesIn.end(); ++edgeIt)
        {
            _nodes.push_back(this->getNode(this->edges[*edgeIt].first));
        }
        for (std::vector<int>::iterator edgeIt = _edgesOut.begin(); edgeIt != _edgesOut.end(); ++edgeIt)
        {
            _nodes.push_back(this->getNode(this->edges[*edgeIt].second));
        }
        return _nodes;
    }

    std::vector<int> Graph::getOutNeighbourKeys(const int& key)
    {
        std::vector<int> _keys;
        std::vector<int> _edges = this->findOutEdges(key);
        for (std::vector<int>::iterator edgeIt = _edges.begin(); edgeIt != _edges.end(); ++edgeIt)
        {
            _keys.push_back(this->edges[*edgeIt].second);
        }
        return _keys;
    }

    std::vector<int> Graph::getInNeighbourKeys(const int& key)
    {
        std::vector<int> _keys;
        std::vector<int> _edges = this->findInEdges(key);
        for (std::vector<int>::iterator edgeIt = _edges.begin(); edgeIt != _edges.end(); ++edgeIt)
        {
            _keys.push_back(this->edges[*edgeIt].first);
        }
        return _keys;
    }

    std::vector<int> Graph::getNeighbourKeys(const int& key)
    {
        std::vector<int> _keys;
        std::vector<int> _edgesIn = this->findInEdges(key);
        std::vector<int> _edgesOut = this->findOutEdges(key);
        for (std::vector<int>::iterator edgeIt = _edgesIn.begin(); edgeIt != _edgesIn.end(); ++edgeIt)
        {
            _keys.push_back(this->edges[*edgeIt].first);
        }
        for (std::vector<int>::iterator edgeIt = _edgesOut.begin(); edgeIt != _edgesOut.end(); ++edgeIt)
        {
            _keys.push_back(this->edges[*edgeIt].second);
        }
        return _keys;
    }