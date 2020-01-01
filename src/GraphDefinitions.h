// save some definitions that are needed for all Graph classes
#include <map>


#ifndef _GRAPHDEFINITIONS_
#define _GRAPHDEFINITIONS_

using namespace std;

// define a struct that contains all the necessary information on the edges
struct Edge {
    double capacity; // what is the capacity of the node
    double flow; // flow at the node; same as derivative
    double tension; // tension on that node; saved as lambda * t
    double lambda; // lambda when the tension was saved
};

struct MaxFlowEdge {
    int to;
    Edge* edgePtr; // pointer to the edge
    Edge* edgePtrBack; // pointer to the backwards edge
};

typedef map<int,Edge*> Node;
typedef map<int,Node> Nodes;

typedef vector<MaxFlowEdge> MaxFlowNode;
typedef vector<MaxFlowNode> MaxFlowNodes;

// const int whiteCode=-3;
const int startCode = -4;
const int source = 0;
const int sink = 1;
const int emptyCode = -5;


#endif
