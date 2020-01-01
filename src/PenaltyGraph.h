/**********************************************************************
***
*** Implements a graph of all the penalties in the fused lasso problems
*** Groups of fused variables will be treated by taking subgraphs
*** 
***********************************************************************/
#ifndef _PENALTYGRAPH_
#define _PENALTYGRAPH_

#include <set>
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include "GeneralFunctions.h"
#include "GraphDefinitions.h"
#include "MaxFlowGraph.h"
#include <R.h>
#include <Rinternals.h>

using namespace std;


class PenaltyGraph
{
public:
    Nodes nodes;
    
    // helper function for subGraph that copies the pointers to the edges in the
    // subgraph, in nodePull the pull on each of the nodes will be saved (needed for the source and sink later)
    void subGraphGetEdges(MaxFlowGraph& m, list<pair<int,double> >& nodePull);
    // helper function for subGraph that generates the source and sink node, nodePull gets deleted
    void subGraphSourceSink(MaxFlowGraph& m, list<pair<int, double> >& nodePull);
    

    // adds an edge to the graph (in both nodes); only intended for use at the start of the 
    // algorithm; will set tension and lambda to 0; flow to the sign (corrected for direction b/c of 2 nodes)
    // and sets the capacity to 1 in the direction of positive sign and infinity in the other (may not be necessary)
    void addEdge(const int from, const int to, const int sign); 
    
    // given a set of nodes, return the set of nodes these nodes are connected to
    // excluding the nodes in the input set itself; these are the nodes a subgraph is connected to
    set<int> connectedTo(const set<int>& subNodes);
    
    // given a list of nodes, the function generates a graph containing a source and a sink
    // the edges in the graphs are pointers to the edges in the PenaltyGraph object, so that
    // operations in the MaxFlowGraph are automatically stored in the PenaltyGraph as well
    MaxFlowGraph* subGraph(const set<int>& subNodes);
    
    // get the sign of the flow between two groups of nodes;
    // this function will only look for the first edge that links the two groups of nodes and
    // take its sign; it will not check that this is consistent over all edges (it will be if 
    // the nodes are from different groups)
    int flowSignBetweenGroups(const set<int>& nodes1, const set<int>& nodes2);
    
    // constructor that uses an R object to build the graph
    // the object is a list; the first element is a vector with the number of the nodes
    // the second element is a list with elements that are vectors of nodenumbers the nodes are
    // connected to
    // startValue is the starting value for the nodes; exact value not important, only ordering
    PenaltyGraph(SEXP connList, SEXP startValue);
    PenaltyGraph(vector<int> nodeNumbers, vector<list<int> > conn, vector<double> startValues);
    PenaltyGraph(){};
    
    // destructor that frees all the edges 
    ~PenaltyGraph();

    // returns the maximum value of the node numbers
    int getMaxNodeNum();
    
    // get a set with all the nodes
    set<int> allNodes();
    
    // prints out the whole graph; used for troubleshooting
    void printGraph(ostream& outStream);
};

#endif
