#include "MaxFlowGraph.h"
#include "GeneralFunctions.h"
#include <queue>
#include <fstream>
#include <sstream>
#include <R.h>
#include <Rinternals.h>


void MaxFlowGraph::addEdgeCap(int from, int to, double capacity)
{
    Edge* e1 = new(Edge);
    Edge* e2 = new(Edge);
    
    // initialize some values
    e1->tension=0;
    e1->lambda=0;
    e2->tension=0;
    e2->lambda=0;
    e1->flow = 0;
    e2->flow = 0;
    
    // set the capacity; capacity in one direction, 0 in the other
    e1->capacity=capacity;
    e2->capacity=0;
    
    // add the elements to the graph
    addEdgeOneWay(from, to, e1, e2);
    addEdgeOneWay(to, from, e2, e1);
}

void MaxFlowGraph::addEdgeOneWay(int from, int to, Edge *edgePtr, Edge *edgePtrBack)
{
    int foo;
    MaxFlowEdge bar;
    foo = nodes[from].size();
    nodes[from].resize(foo+1);
    bar.to = to;
    bar.edgePtr = edgePtr;
    bar.edgePtrBack = edgePtrBack;
    nodes[from][foo]=bar;
}

void MaxFlowGraph::deleteAllEdges(int nodeNum)
{
    MaxFlowNode::iterator edgeIt;
    
    // go through all the nodes in the node nodeNum and delete them
    for(edgeIt = nodes[nodeNum].begin(); edgeIt!=nodes[nodeNum].end(); ++edgeIt)
    {
        delete(edgeIt->edgePtr);
        delete(edgeIt->edgePtrBack);
    }
    nodes[nodeNum].clear();
}



pair<int,int> MaxFlowGraph::addSpecialSourceSink(const vector<double>& overFlow)
{
    int newSource = nodes.size();
    int newSink = nodes.size()+1;
    nodes.resize(nodes.size()+2);
    for(unsigned int i=0; i< overFlow.size(); ++i)
    {
        if(overFlow[i]>0)
        {
            addEdgeCap(newSource, i, overFlow[i]);
        }
        else if(overFlow[i]<0)
        {
            addEdgeCap(i,newSink, -overFlow[i]);
        }
    }
    return(pair<int,int>(newSource, newSink));
}

void MaxFlowGraph::removeSpecialSourceSink(const vector<double>& overFlow, const int newSource, const int newSink)
{
    int numEdges;
    // first go through all other nodes and delete the appropriate edge, which by construction has to be the last one
    for(unsigned int i=0; i< overFlow.size(); ++i)
    {
        if(overFlow[i]!=0)
        {
            numEdges = nodes[i].size();
            nodes[i].erase(nodes[i].begin()+(numEdges-1));
        }
    }
    // now delete the appropriate edges
    deleteAllEdges(newSource);
    deleteAllEdges(newSink);
    // now delete the nodes itself; do it with max and min because numbering changes after erasing
    nodes.erase(nodes.begin()+Max(newSource,newSink));
    nodes.erase(nodes.begin()+Min(newSource,newSink));
}













vector<int> MaxFlowGraph::distance(int start, bool from) 
{
    vector<int> dist(nodes.size(), nodes.size()); // initialize the parent vector
    queue<int> next;
    int u; // number of node that is currently being looked at
    Edge* e;
    MaxFlowNode::iterator edgeIt;
    
    dist[start]=0; 
    next.push(start); // start is the starting point
    
    while(!next.empty())
    {
        u=next.front();
        next.pop();

        // checking node u
        // go through all edges in node u;
        for(edgeIt= nodes[u].begin(); edgeIt!=nodes[u].end(); ++edgeIt)
        {
            if(from) // distance from the last node; depends which direction to check capacity
            {
                e=edgeIt->edgePtr;
            }
            else // to the last node
            {
                e = edgeIt->edgePtrBack; // this way the distance when flowing towards start is calculated
            }
            if(e->flow < e->capacity-tolerance) // is reachable
            {
                if(dist[edgeIt->to]>dist[u]+1) // distance can be updated
                {
                    dist[edgeIt->to]=dist[u]+1;
                    next.push(edgeIt->to);
                }
            }
        }
    }
    return(dist);
}

bool MaxFlowGraph::push(int from, MaxFlowEdge &e, const int sourceNode, const int sinkNode)
{
    bool isToActive;
    double foo = Min(exFlow[from], e.edgePtr->capacity - e.edgePtr->flow);
    e.edgePtr->flow +=foo;
    e.edgePtrBack->flow -= foo;
    exFlow[from]-=foo;
    // save if to is not already an active node
    isToActive = (exFlow[e.to]>tolerance);
    exFlow[e.to]+=foo;
    if(!isToActive && e.to!=sourceNode && e.to!=sinkNode) // add it to the active nodes if not active already and not source or sink
    {
        insertActiveNode(e.to);
    }
    return(exFlow[from]>tolerance); // return true if the node is still active
}

int MaxFlowGraph::findDist(int nodeNum)
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    int newDist = infiniteInt;
    
    nodeIt = nodes.begin()+nodeNum;
    for(edgeIt = nodeIt->begin(); edgeIt!=nodeIt->end(); ++edgeIt)
    {
        // check that there is residual capacity left
        if(edgeIt->edgePtr->flow < edgeIt->edgePtr->capacity - tolerance)
        {
            newDist = Min(newDist, dist[edgeIt->to]+1);
        }
    }
    return(newDist);
}


void MaxFlowGraph::preprocess(const int sourceNode, const int sinkNode)
{
    // compute the distance label from the sink using the current flow
    dist = distance(sinkNode);
    
    // clean the active nodes;
    activeByDist.assign(2*nodes.size()+1, list<int>(0));
    level=-1;
    
    // set all the excess flows to 0
    exFlow.assign(nodes.size(),0);
    
    // go through all nodes connected to the source and push the maximum flow
    MaxFlowNode::iterator edgeIt;
    for(edgeIt = nodes[sourceNode].begin(); edgeIt!=nodes[sourceNode].end(); ++edgeIt)
    {
        // what is the excess in the connected node
        exFlow[edgeIt->to] = edgeIt->edgePtr->capacity - edgeIt->edgePtr->flow;
        exFlow[sourceNode] -=exFlow[edgeIt->to];
        // set flow to maximum
        edgeIt->edgePtr->flow = edgeIt->edgePtr->capacity;
        edgeIt->edgePtrBack->flow = -edgeIt->edgePtr->capacity;
        // if positive excess flow, insert it as an active node
        if(exFlow[edgeIt->to]>tolerance)
        {
            insertActiveNode(edgeIt->to);
        }
    }
    dist[sourceNode] = nodes.size();
}

bool MaxFlowGraph::pushRelabel(const int i, const int sourceNode, const int sinkNode)
{
    // go through all edges from node i and look for admissible ones
    // that is d(i)= d(j)+1 and flow(i,j) < cap(i,j)
    bool didPush = false; // saves if a successfull push has occured
    MaxFlowNode::iterator edgeIt;
    for(edgeIt=nodes[i].begin(); edgeIt!=nodes[i].end(); ++edgeIt)
    {
        // check if the edge is admissible
        if((dist[i]==dist[edgeIt->to]+1) && (edgeIt->edgePtr->capacity > edgeIt->edgePtr->flow + tolerance))
        {
            didPush=true;
            if(!push(i,*edgeIt, sourceNode, sinkNode)) // node became inactive
            {
                return(false);
            }
        }
    }
    if(!didPush) // no successfull push was possible; relabel
    {
        dist[i] = findDist(i);
    }
    
    return(true); // node still active as there is still an excess
}


bool MaxFlowGraph::getLargestActiveNode(int& nodeNum)
{
    // check that there is an active node with distance level
    if(level<0) // nothing in here
    {
        return(false);
    }
    if(activeByDist[level].empty()) // no there isn't, step down until there is
    {
        do
        {
            --level;
        }
        while((level>=0) && (activeByDist[level].empty()) ); // check first if level still >=0, then if list is empty
        if(level<0) // no active nodes left
        {
            return(false);
        }
    }
    // save the largest element and delete it;
    nodeNum = activeByDist[level].front();
    activeByDist[level].pop_front();
    return(true); // saved active node
}

void MaxFlowGraph::insertActiveNode(const int nodeNum)
{
    // check if new node is the largest
    if(dist[nodeNum]>level)
    {
        level=dist[nodeNum];
    }
    // insert the node
    activeByDist[dist[nodeNum]].push_back(nodeNum);
}


bool MaxFlowGraph::checkSourceMaxedOut(const int sourceNode)
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    nodeIt = nodes.begin()+sourceNode;
    for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt)
    {
        if(edgeIt->edgePtr->flow < edgeIt->edgePtr->capacity-tolerance)
        {
            return(false); // not maxed out
        }
    }
    
    return(true); // all edges are maxed out
}



double MaxFlowGraph::maxFlowFromSource(int sourceNode)
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    double foo=0;
    
    nodeIt = nodes.begin()+sourceNode;
    for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt)
    {
        foo+=edgeIt->edgePtr->capacity;
    }
    return(foo);
}

double MaxFlowGraph::currentFlowFromSource(int sourceNode)
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    double foo=0;
    
    nodeIt = nodes.begin()+sourceNode;
    for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt)
    {
        foo+=edgeIt->edgePtr->flow;
    }
    return(foo);
}



void MaxFlowGraph::setCapacity()
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    for(nodeIt = nodes.begin()+2; nodeIt != nodes.end(); ++nodeIt) // iterate through all nodes, but jump over source and sink
    {
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) // iterate through all edges
        {
            // check that it is not the source or sink
            if((edgeIt->to!=source) && (edgeIt->to!=sink))
            {
                // set capacity only to 1 if maximal tension reached; otherwise infinity
                if((edgeIt->to!=source) && (edgeIt->to!=sink))
                {
                // set capacity only to 1 if maximal tension reached; otherwise infinity
                    if(RelDif(edgeIt->edgePtr->tension, edgeIt->edgePtr->lambda)> tolerance)
                    {
                        edgeIt->edgePtr->capacity = infinite;
                    }
                    else
                    {
                        edgeIt->edgePtr->capacity = 1;
                    }
                }
            }
        }
    }
}



void MaxFlowGraph::setCapacityTo1()
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    for(nodeIt = nodes.begin()+2; nodeIt != nodes.end(); ++nodeIt) // iterate through all nodes, but jump over source and sink
    {
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) // iterate through all edges 
        {
            if((edgeIt->to!=source) && (edgeIt->to!=sink))
            {
                edgeIt->edgePtr->capacity = 1;
            }
        }
    }
}


void MaxFlowGraph::setCapacityProportional(double factor)
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    Edge* e;
    
    for(nodeIt = nodes.begin()+2; nodeIt != nodes.end(); ++nodeIt) // iterate through all nodes, but jump over source and sink
    {
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) // iterate through all edges
        {
            // check that it is not the source or sink
            if((edgeIt->to!=source) && (edgeIt->to!=sink))
            {
                // set capacity only to 1 if maximal tension reached; otherwise infinity
                e=edgeIt->edgePtr;
                if(RelDif(e->tension, e->lambda)> tolerance)
                {
                    e->capacity = 1+ factor*RelDif(e->lambda,e->tension);
                }
                else
                {
                    e->capacity = 1;
                }
            }
        }
    }
}


void MaxFlowGraph::updateCapacity(const double newLambda, vector<double>& overFlow)
{
    // initialize the vector with 0s
    overFlow.assign(nodes.size(),0);
    
    // adjust the edges if necessary
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    int nodeNum;
    for(nodeIt = nodes.begin()+2, nodeNum=2; nodeIt != nodes.end(); ++nodeIt, ++nodeNum) // iterate through all nodes, but jump over source and sink
    {
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) // iterate through all edges 
        {
            if((edgeIt->to!=source) && (edgeIt->to!=sink))
            {
                if(edgeIt->edgePtr->capacity==1)
                {
                    if(edgeIt->edgePtr->tension < edgeIt->edgePtr->lambda -tolerance) // no flow adjustments necessary
                    {
                        edgeIt->edgePtr->capacity=infinite;
                    }
                }
                else if(edgeIt->edgePtr->capacity>1)
                {
                    if(edgeIt->edgePtr->tension >= edgeIt->edgePtr->lambda-tolerance) // adjust the flow and capacity downwards
                    {
                        edgeIt->edgePtr->capacity =1;
                        if(edgeIt->edgePtr->flow >1) // correct the flow downward if necessary
                        {
                            overFlow[nodeNum] +=edgeIt->edgePtr->flow -1;
                            overFlow[edgeIt->to] -= edgeIt->edgePtr->flow -1;
                            edgeIt->edgePtr->flow=1;
                            edgeIt->edgePtrBack->flow=-1;
                        }
                    }
                }
            }
        }
    }
}



void MaxFlowGraph::setFlowTo0()
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    
    for(nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt) // iterate through all nodes
    {
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) // iterate through all edges
        {
            edgeIt->edgePtr->flow = 0;
        }
    }
}








MaxFlowGraph::MaxFlowGraph(const set<int>& graphNodes):nodes(graphNodes.size()+2,MaxFlowNode(0)), exFlow(graphNodes.size()+2,0), dist(graphNodes.size()+2, graphNodes.size()+2), activeByDist(2*graphNodes.size()+9, list<int>(0)), nodeMapIntToExt(graphNodes.size()+2,-1)
{
    int curIntNode=2;
    set<int>::iterator setIt;
    //set the mapping of the internal to the external nodes
    for(setIt = graphNodes.begin(); setIt!=graphNodes.end(); ++setIt)
    {
        nodeMapIntToExt[curIntNode]=*setIt;
        nodeMapExtToInt[*setIt]=curIntNode;
        ++curIntNode;
    }
    groupDeriv=0;
    lambda=0;
    
}




// Do the preflow-push algorithm to find the maximal flow, using the currently saved flow;
// initialization to 0 has to be done seperately by hand
bool MaxFlowGraph::findMaxFlow(const int sourceNode, const int sinkNode)
{
    int activeNodeNum;
    preprocess(sourceNode, sinkNode);
    while(getLargestActiveNode(activeNodeNum))
    {
        if(pushRelabel(activeNodeNum,sourceNode, sinkNode)) // node remains active; if inactive, does not need to be replaced
        {
            insertActiveNode(activeNodeNum);
        }
    }

    return(checkSourceMaxedOut(sourceNode));
}


double MaxFlowGraph::validUntil()
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    Edge *e;

    double validLambda=infinite;
    double foo, offset;
    
    // go through all edges except those connected to the source or the sink
    for(nodeIt = nodes.begin()+2; nodeIt != nodes.end(); ++nodeIt)
    {
        for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt)
        {
            // check that it is not pointing to the source or sink
            if((edgeIt->to!=source) && (edgeIt->to!=sink))
            {
                if(edgeIt->edgePtr->flow > 1+tolerance) // calculate when tension will hit lambda (which is also changing)
                {
                    e = edgeIt->edgePtr;
                    {
                        offset = ((e->lambda - e->tension)/(e->flow-1));
                        if(offset < 0) // should not happen
                        {
                            e->tension = e->lambda;
                            edgeIt->edgePtrBack->tension = -e->lambda;
                        }
                        else
                        {
                            foo = e->lambda + offset;
                            validLambda = Min(validLambda, foo);
                        }
                    }
               }
            }
        }
    }
    if(validLambda==infinite)
    {
        validLambda=neverSplit;
    }
    return(validLambda);
}





















void MaxFlowGraph::updateTension(const double newLambda)
{
    MaxFlowNodes::iterator nodeIt;
    MaxFlowNode::iterator edgeIt;
    Edge* e;
    
    {
        for(nodeIt = nodes.begin()+2; nodeIt != nodes.end(); ++nodeIt)
        {
            for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt)
            {
                // check that it is not pointing to the source or sink
                if((edgeIt->to!=source) && (edgeIt->to!=sink))
                {
                    e = edgeIt->edgePtr;
                    e->tension += e->flow*(newLambda-e->lambda);
                    e->lambda = newLambda;
                }
            }
        }
    }
    lambda = newLambda;
}



double MaxFlowGraph::calcTensionChange(const double lambda)
{
    // first update the tensions
    updateTension(lambda);
    
    // set the flow to 0 as preparation for calculating tension derivatives
    setFlowTo0();
    
    // first set the capacity to 1, find an optimal flow to see if the group will ever split
    setCapacityTo1();
    if(findMaxFlow(source, sink)) // if returns true, then all source edges are at capacity
    {
        return(neverSplit);
    }
    else // increase capacities to maximum and maximize flow again
    {
//            cout << " Full Capacity" << endl;
        setCapacity();
        if(findMaxFlow(source, sink))
        {
            return(validUntil());
        }
        else
        {
            return(splitNow);
        }
    }
}



double MaxFlowGraph::calcTensionChangeUpdate(const double lambda)
{
    stringstream updateOutput;
    // update the tension in the graph
    updateTension(lambda);
    
    // update the capacity, tension update called above
    vector<double> overFlow;
    updateCapacity(lambda, overFlow);
    
    // build the new source and sink node
    pair<int,int> newNodes = addSpecialSourceSink(overFlow);
    
    // do the maximum flow
    bool result = findMaxFlow(newNodes.first, newNodes.second);
    // remove the new source and sink node
    removeSpecialSourceSink(overFlow, newNodes.first, newNodes.second);
    
    // if found a new valid flow, return how long it lasts
    if(result)
    {
        return(validUntil());
    }
    else // there is no valid maximum flow anymore; call a regular maxFlow from 0 to find out where to split
    {
        setFlowTo0();
//        ofstream outStream("TestOutput3.txt", ios::out);
        result = findMaxFlow(source, sink);
        if(result) // this should not happen, stop with an error (can be thrown out later because it should be unncecessary
        {
/*
            cout << "Found a possible flow in calcTensionChangeUpdate although deemed impossible" << endl;
            char foo[256];
            while(!updateOutput.eof())
            {
                updateOutput.getline(foo,256);
                cout << foo << endl;
            }
            cout << "Graph after new flow" << endl << foo << endl;

            printGraph(cout);
            exit(1);
*/
        }
        return(splitNow);
    }
}


double MaxFlowGraph::calcTensionChangeProportional(const double lambda)
{
    // first update the tensions
    double currentFlow, currentFactor, deltaFactor, deltaFlow, maxFlow;
    updateTension(lambda);
    
    // set the flow to 0 as preparation for calculating tension derivatives
    setFlowTo0();
    
    // first set the capacity to 1, find an optimal flow to see if the group will ever split
    setCapacityTo1();
    if(findMaxFlow(source, sink)) // if returns true, then all source edges are at capacity
    {
        return(neverSplit);
    }
    
    currentFlow=currentFlowFromSource(source);
    maxFlow = maxFlowFromSource(source);
    currentFactor=(maxFlow-currentFlow)/currentFlow/2;// a factor of 1 raises the capacity at most by 1;
    deltaFactor = currentFactor;
    deltaFlow = currentFlow;
    setCapacityProportional(currentFactor);
    while(!findMaxFlow(source, sink))
    {
        deltaFlow = currentFlowFromSource(source) - currentFlow;
        currentFlow += deltaFlow;
        // now calculate the new necessary increase in the factor
        deltaFactor = deltaFactor * (maxFlow-currentFlow)/deltaFlow;
        currentFactor+=deltaFactor;
        
        if(deltaFlow < tolerance){ // no increase anymore; won't finish
            return(splitNow);}
        
        setCapacityProportional(currentFactor);
    }
    return(validUntil());
}



set<int> MaxFlowGraph::reachableFromSource(const int sourceNode)
{
    set<int> reachable;
    vector<int> distReach=distance(sourceNode,true); // get the distance from the source; nodes.size() means that it can't be
    // reached

    // reachable nodes are all nodes in parent except source and sink
    for(unsigned int i=2; i !=distReach.size(); ++i)
    {
      if(((unsigned int) distReach[i])<nodes.size())
        {
            // node is reachable; insert node in external notation
            reachable.insert(nodeMapIntToExt[i]); 
        }
    }
    return(reachable);
}


set<int> MaxFlowGraph::getComplement(const set<int> &x)
{
    set<int> complement;
    map<int,int>::const_iterator nodeIt;
    for(nodeIt=nodeMapExtToInt.begin(); nodeIt!=nodeMapExtToInt.end(); ++nodeIt)
    {
        // do not have to check for source and sink nodes as they are not included in nodeMap variables
        if((x.count(nodeIt->first)==0)) // not in x, so add to complement
        {
            complement.insert(nodeIt->first);
        }
    }
    return(complement);
}


set<int> MaxFlowGraph::allNodes()
{
    set<int> all;
    map<int,int>::const_iterator nodeIt;
    for(nodeIt=nodeMapExtToInt.begin(); nodeIt!=nodeMapExtToInt.end(); ++nodeIt)
    {
        // do not have to check for source and sink nodes as they are not included in nodeMap variables
        all.insert(nodeIt->first);
    }
    return(all);
}




void MaxFlowGraph::clear()
{
    // delete all edges go to and from the source and sink nodes
    deleteAllEdges(source);
    deleteAllEdges(sink);

    // delete all the pointers to the other edges
    nodes.clear();
    groupDeriv = 0;
}

MaxFlowGraph::~MaxFlowGraph()
{
    // delete all edges go to and from the source and sink nodes
    deleteAllEdges(source);
    deleteAllEdges(sink);

    nodeMapExtToInt.clear();
    nodeMapIntToExt.clear();
    nodes.clear();
}
    
    
/*
void MaxFlowGraph::printGraph(ostream& outStream)
{
    int edgeNum, nodeNum;
    MaxFlowEdge me;
    
    outStream << "Group movement: " << groupDeriv << endl;
    for(nodeNum = 0; nodeNum != nodes.size(); ++nodeNum)
    {
        // right name for the node; internal to external mapping with special case for source and sink
        if(nodeNum==source)
        {
            outStream << "Node Number: Source " << nodeNum << endl;
        }
        else if(nodeNum==sink)
        {
            outStream << "Node Number: Sink " << nodeNum << endl;
        }
        else
        {
            outStream << "Node Number: " << nodeNum << ", " << nodeMapIntToExt[nodeNum] << endl;
        }
        // write the excess flow and the distance
        outStream << "Excess Flow: "<< exFlow[nodeNum] << " Distance: " << dist[nodeNum] << endl;
        
        outStream << "Edges:" << endl;
        for(edgeNum = 0; edgeNum != nodes[nodeNum].size(); ++edgeNum)
        {
            me = nodes[nodeNum][edgeNum];
            if(me.to==source)
            {
                outStream << "To: Source";
            }
            else if(me.to==sink)
            {
                outStream << "To: Sink";
            }
            else
            {
                outStream << "To: " << me.to;
            }
            double bar = me.edgePtr->lambda + ((me.edgePtr->lambda - me.edgePtr->tension)/(me.edgePtr->flow-1));
            outStream.precision(15);
            outStream << " Cap: " << me.edgePtr->capacity <<
                " Flow: " << me.edgePtr->flow << " Tension: " << me.edgePtr->tension << " Lambda: " << 
                 me.edgePtr->lambda;
#ifdef _DEBUG2_
            outStream << "Diff: " << me.edgePtr->tension - me.edgePtr->lambda << 
                 "Delta:" << ((me.edgePtr->lambda - me.edgePtr->tension)/(me.edgePtr->flow-1))<< "Valid: " << bar <<  "Diff:" << bar - me.edgePtr->lambda << "Eligible: " << (me.edgePtr->flow > 1);
#endif
            outStream.precision(5);
            outStream << endl;
        }
        outStream << endl;
    }
    
//    map<int,int>::iterator MI;
//    for(MI = nodeMapExtToInt.begin(); MI!=nodeMapExtToInt.end(); ++MI)
//    {
//        outStream << MI->first << " " << MI->second << endl;
//    }
    
    outStream << endl;
}


void MaxFlowGraph::printActiveNodes(ostream& outStream)
{
    int i;
    list<int>::iterator listIt;
    
    for(i=0; i<activeByDist.size(); ++i)
    {
        for(listIt = activeByDist[i].begin(); listIt !=activeByDist[i].end(); ++listIt)
        {
            outStream << "Dist: " << i << " Node: " << *listIt << endl;
        }
    }
}
*/


/*
#include <iostream>
int main(int argc, char** argv)
{
    set<int> foo;
    foo.insert(2);
    foo.insert(3);
    foo.insert(4);
    foo.insert(5);
    foo.insert(6);
    foo.insert(7);
    MaxFlowGraph m(foo);
    
    // insert the graph
    
    
    m.addEdgeCap(3,2,1);
    m.addEdgeCap(3,4,infinite);
    
    m.addEdgeCap(4,5,infinite);
    m.addEdgeCap(6,5,infinite);
    m.addEdgeCap(5,7, infinite);
    
    
    m.addEdgeCap(7,sink,4.0/3);
    m.addEdgeCap(3,sink,1.0/3);
    m.addEdgeCap(2,sink,4.0/3);
    m.addEdgeCap(source,6,2.0/3);
    m.addEdgeCap(source,5,2.0/3);
    m.addEdgeCap(source,4,5.0/3);
    
    // adjust a few capacities
    m.nodes[2][1].edgePtr->capacity = infinite;
    m.nodes[4][1].edgePtr->capacity = infinite;
    m.nodes[5][1].edgePtr->capacity = infinite;
    m.nodes[5][3].edgePtr->capacity = infinite;
    m.nodes[7][0].edgePtr->capacity = infinite;
    
    
    m.printGraph(cout);
    
    m.setFlowTo0();
    m.findMaxFlow(source, sink);
    m.printGraph(cout);
    
    
    return(0);
}
*/
/*
    m.addEdgeCap(source,2,5);
    m.addEdgeCap(source,3,2);
    m.addEdgeCap(source,4,1);
    
    m.addEdgeCap(2,5,5);
    m.addEdgeCap(2,6,5);
    m.addEdgeCap(2,7,5);
    
    m.addEdgeCap(3,5,6);
    m.addEdgeCap(3,6,6);
    m.addEdgeCap(3,7,6);
    
    m.addEdgeCap(4,5,12);
    
    m.addEdgeCap(5,sink,10);
    m.addEdgeCap(6,sink,15);
    m.addEdgeCap(7,sink,15);
*/
