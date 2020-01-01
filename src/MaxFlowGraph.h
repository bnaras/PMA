#ifndef _MAXFLOWGRAPH_
#define _MAXFLOWGRAPH_

//#include <iostream>
#include <utility>
#include <list>
#include <set>
#include <vector>
#include "GraphDefinitions.h"
#include "GeneralFunctions.h"

const double neverSplit=-1;
const double splitNow=-2;



class MaxFlowGraph
{
    MaxFlowNodes nodes;
    // variable for the preflow-push algorithm
    vector<double> exFlow; // saves the excess flow
    vector<int> dist; // distance label from the sink
    vector<list<int> > activeByDist; // active nodes by distance
    int level; // largest distance of an active node
    
    // the following variables map internal to external nodes; source and sink always excluded
    map<int,int> nodeMapExtToInt; // maps the external node numbers to the internal node numbers
    vector<int> nodeMapIntToExt; // maps the internal nodes to the external ones
    double groupDeriv; // mean value of the pull on each node
    double lambda;
    
    friend class PenaltyGraph;

    // The following functions are used for generating a valid flowGraph
    // add an edge to the graph; only the capacity is set in one direction, everything else to 0;
    // usually only used for the source and sink node; from and to are in internal numeration
    void addEdgeCap(int from, int to, double capacity);
    // adds an edge only to the from node by increasing the vector size; notation in internal numeration
    void addEdgeOneWay(int from, int to, Edge *edgePtr, Edge *edgePtrBack);
    // frees the memory for all edges in a node; should only be used on source or sink; internal numeration
    // helper function for the destructor
    void deleteAllEdges(int nodeNum);
    
    // the following functions add and remove special source and sink nodes after all other nodes
    // have been introduced; the default source and sink nodes are treated as regular nodes
    // after calling add, remove has to be called at some point to clean up
    pair<int,int> addSpecialSourceSink(const vector<double>& overFlow); // returns the new source and sink numbers
    void removeSpecialSourceSink(const vector<double>& overFlow, const int newSource, const int newSink);
    
    
    // The following functions are helper functions for finding maximum flows
    
    // function performs a breadth first search; return value is the parent of each node
    // if from is true, then it will calculate the distance from start, otherwise the distance to start
    vector<int> distance(int start, bool from=false);

    // new helper function for preflow-push algorithm needed; returns true if the node is still active after the push
    bool push(int from, MaxFlowEdge& e, const int sourceNode, const int sinkNode);
    // finds the new label for node nodeNum
    int findDist(int nodeNum);
    // preprocess the graph; sets the distances from sink and pushes maximum flow from the source
    void preprocess(const int sourceNode, const int sinkNode);
    // push/relabel procedure; return value indicates if node is still active
    bool pushRelabel(const int i, const int sourceNode, const int sinkNode);
    // return the next largest active node; will return false if there are no more; node will be saved in nodeNum
    bool getLargestActiveNode(int& nodeNum);
    // insert an active node into the data structure
    void insertActiveNode(const int nodeNum);

    // checks if all edges of the source are at maximum (from node point of view); only sensible for the source
    bool checkSourceMaxedOut(const int sourceNode); 

    // gives the maximum amount possible that can come out of the source
    double maxFlowFromSource(int sourceNode);
    // gives the amount currently flowing from the source
    double currentFlowFromSource(int sourceNode);
    
    // functions to initialize capacity, flow and update tension; only work for source=0 and sink=1;
    // not needed for general source and sink node, only for source=0, sink=1
    void setCapacity(); // sets the capacity, excpet for edges to and from source or sink 
    void setCapacityTo1(); // sets the capacity to 1 everywhere except for edges coming or going from the sink or source
    // sets capacity to 1+ proportion to leftover room from tension to lambda
    void setCapacityProportional(double factor);
    
    // change the capacity only at edges where necessary, change the flow to fit the capacity and return 
    // the changes in the overFlow vector given as input
    void updateCapacity(const double newLambda, vector<double>& overFlow);
    
    void setFlowTo0(); // flow set to 0 everywhere
    
    // returns true if all edges from the sourceNode are at full capacity
    bool findMaxFlow(const int sourceNode, const int sinkNode); 
    
    // calculate until when the current flows are valid
    double validUntil();

    // deletes all nodes in the graph; frees the edges from the source and sink
    void clear();

public:
    // constructor for the class; need to give a list with nodes that will be in the graph
    MaxFlowGraph(const set<int>& graphNodes);
    
    // updates the tension values to a higher value of lambda
    void updateTension(const double newLambda); 
    
    // calculates the derivative of the tension; return value gives the lambda value when it has to be recalculated
    // negative values can indicate that it has never to be recalculated or that the group has to be split (see constants above)
    double calcTensionChange(const double lambda);
    // while the previous function can be called any time, this function can only be called after the
    // previous function has been called at least once. However it is particularly efficient
    // at doing flow updates when a tension update is necessary because a boundary has been hit
    double calcTensionChangeUpdate(const double lambda);
    // uses the proportional version of capacity setting to minimize the number of necessary tension calls
    double calcTensionChangeProportional(const double lambda);
    
    // returns a set of nodes that are reachable from the source with current flow
    set<int> reachableFromSource(const int sourceNode=source);
    // get the complement of the set
    set<int> getComplement(const set<int>& x);
    // returns all nodes contained in the flowGraph (except source and sink)
    set<int> allNodes();
        
    // returns the derivative for the whole group;
    // is the derivative of the group mu
    inline double getGroupDeriv() {return(groupDeriv);};
    inline int size() {return(nodes.size()-2);};

    // destruktor for MaxFlowGraph object
    ~MaxFlowGraph();

    // prints out the whole graph; used for troubleshooting
    //void printGraph(ostream& outStream);
    // prints the currently active nodes
    //void printActiveNodes(ostream& outStream);
    
};

#endif
