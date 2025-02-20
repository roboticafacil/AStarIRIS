#ifndef NODE_H
#define NODE_H
class Node
{
public:
    Node();
    virtual ~Node();
    virtual void * getNodeData()=0;
};
#endif