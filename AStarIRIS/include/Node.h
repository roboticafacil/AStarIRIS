#ifndef NODE
#define NODE
class Node
{
public:
    Node();
    virtual ~Node();
    virtual void * getNodeData()=0;
};
#endif