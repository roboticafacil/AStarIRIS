#ifndef EDGE_H
#define EDGE_H
class Edge
{
public:
    Edge();
    virtual ~Edge();
    virtual void* getEdgeData() = 0;
};
#endif