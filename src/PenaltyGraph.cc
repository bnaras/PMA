#include "PenaltyGraph.h"
#include <vector>


void PenaltyGraph::addEdge(int from, int to, int sign)
{
    Edge* e1 = new(Edge);
    Edge* e2 = new(Edge);
    
    // initialize some values
    e1->tension=0;
    e1->lambda=0;
    e2->tension=0;
    e2->lambda=0;
    
    // set the flow (also works as derivative)
    e1->flow = sign;
    e2->flow = -sign;
    
    // set the capacity
    if(sign==1)
    {
        e1->capacity=sign;
        e2->capacity=infinite;
    }
    else if(sign==-1)
    {
        e1->capacity=infinite;
        e2->capacity=-sign;
    }
    else if(sign==0)
    {
         e1->capacity = 1;
         e2->capacity= 1;
    }
    else
    {
        throw("Wrong sign given in addEdge");
    }
    
    nodes[from][to]=e1;
    nodes[to][from]=e2;
}



set<int> PenaltyGraph::connectedTo(const set<int>& subNodes)
{
    set<int> conn; // saves the nodes it is connected to
    
    set<int>::const_iterator setIter;
    Nodes::iterator nodeIt;
    Node::iterator edgeIt;
    
    //step through all nodes in the set
    for(setIter=subNodes.begin(); setIter!=subNodes.end(); ++setIter)
    {
        nodeIt = nodes.find(*setIter);
        if(nodeIt!=nodes.end()) // node exists
        {
            // walk through the edges
            for(edgeIt = nodeIt->second.begin(); edgeIt!= nodeIt->second.end(); ++edgeIt) 
            {
                if(subNodes.count(edgeIt->first)==0) // edge to an outside node
                {
                    conn.insert(edgeIt->first); // insert it into the set of connected nodes
                }
            }
        }
    }
    return(conn);
}






void PenaltyGraph::subGraphGetEdges(MaxFlowGraph& m, list<pair<int,double> >& nodePull)
{
    map<int,int>::iterator MI; // iterate through the edges in the subgraph
    Node::iterator edgeIt; // iterate through the edges in each node; check for each if in subgraph
    Nodes::iterator nodeIt; // save the position of the node the algorithm is working on
    pair<int,double> foo; // used to save intermediate results
    int fromNodeIntNum, toNodeIntNum; // saves the internal number of the current node
    Edge *edgePtr, *edgePtrBack;
    
    // go through all the nodes in the subgraph
    for(MI=m.nodeMapExtToInt.begin(); MI!=m.nodeMapExtToInt.end(); ++MI)
    {
        // get iterator for the node in PenaltyGraph
        nodeIt = nodes.find(MI->first);
        // foo is used to calculate the pull on the current node; stored in internal notation
        foo.first = MI->second;
        foo.second=0;
        fromNodeIntNum = MI->second;
        // cycle through all the edges in the current node
        for(edgeIt=nodeIt->second.begin(); edgeIt!=nodeIt->second.end(); ++edgeIt)
        {
            // for each edge check if it is in the subgraph (targetnode of the edge)
            if(m.nodeMapExtToInt.count(edgeIt->first))
            {
                // only check edges for which the target node has higher number than origin node
                // done, as edge both ways is included immediately
                if(edgeIt->first > nodeIt->first)
                {
                    // copy the pointer to the edge
                    toNodeIntNum = m.nodeMapExtToInt[edgeIt->first];
                    edgePtr = edgeIt->second;
                    edgePtrBack = nodes[edgeIt->first][nodeIt->first];
                    // add it one way and the other way
                    m.addEdgeOneWay(fromNodeIntNum, toNodeIntNum,edgePtr,edgePtrBack);
                    m.addEdgeOneWay(toNodeIntNum, fromNodeIntNum, edgePtrBack, edgePtr);
                }
            }
            else // count the pull on the node
            {
                foo.second-= edgeIt->second->flow;
            }
        }
        m.groupDeriv+=foo.second;
        nodePull.push_front(foo); // total pull on the node
    }
    m.groupDeriv/=m.nodeMapExtToInt.size(); // calculate the derivate for the group value as the mean of the pulls
}

void PenaltyGraph::subGraphSourceSink(MaxFlowGraph& m, list<pair<int, double> >& nodePull)
{
    pair<int,double> foo; // used to save intermediate results
    double netPull; // helper variable
    
    // now generate the source and sink node in MaxFlowGraph with the appropriate edges
    while(!nodePull.empty())
    {
        // get the information about the first node
        foo = nodePull.front();
        nodePull.pop_front();
        
        // add the edge for the source or sink node
        netPull = foo.second - m.groupDeriv;
        if(netPull>0)
        {
            m.addEdgeCap(source, foo.first, netPull);
        }
        else if(netPull<0)
        {
            m.addEdgeCap(foo.first, sink, -netPull);
        }
    }
}

MaxFlowGraph* PenaltyGraph::subGraph(const set<int>& subNodes)
{
    MaxFlowGraph* m = new MaxFlowGraph(subNodes);
    list<pair<int,double> >  nodePull; // saves the pull on the nodes; stored with nodeNumbers in internal notation

    subGraphGetEdges(*m, nodePull); // remember, nodePull will be changed, so will m
    subGraphSourceSink(*m, nodePull);

    return(m);
}


int PenaltyGraph::flowSignBetweenGroups(const set<int>& nodes1, const set<int>& nodes2)
{
    set<int>::const_iterator setIter;
    Nodes::iterator nodeIt;
    Node::iterator edgeIt;
    // go through all the nodes in nodes1, check if the edge leads to nodes2, if yes, return the sign
    // of the flow
    for(setIter=nodes1.begin(); setIter!=nodes1.end(); ++setIter)
    {
        nodeIt = nodes.find(*setIter);
        for(edgeIt = nodeIt->second.begin(); edgeIt!=nodeIt->second.end(); ++edgeIt)
        {
            // check if the edge leads to a node in nodes2
            if(nodes2.count(edgeIt->first)) // yes
            {
                return(signum(edgeIt->second->flow));
            }
        }
    }
    // this should not happen
    throw("Asked for sign of flow between groups of nodes that are unconnected in flowSignBetweenGroups of PenaltyGraph");
}



void PenaltyGraph::printGraph(ostream& outStream)
{
    Node::iterator edgeIt; // iterator over the edges
    Nodes::iterator nodeIt; // iterator over the nodes
    
    for(nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
    {
        outStream << "Node Number: " << nodeIt->first << endl;
        outStream << "Edges:" << endl;
        for(edgeIt = nodeIt->second.begin(); edgeIt != nodeIt->second.end(); ++edgeIt)
        {
            outStream << "To: " << edgeIt->first << " Cap: " << edgeIt->second->capacity <<
                " Flow: " << edgeIt->second->flow << " Tension: " << edgeIt->second->tension << " Lambda: " << 
                 edgeIt->second->lambda << endl;
        }
        outStream << endl;
    }
    outStream << endl;
}

PenaltyGraph::PenaltyGraph(SEXP connList, SEXP startValue)
{
    SEXP nodeNumbersR = VECTOR_ELT(connList,0); // numbers of nodes (R VECTOR)
    SEXP connR = VECTOR_ELT(connList,1); // connections of each node (R LIST)
    SEXP connOneNode;
    int numberOfNodes = LENGTH(nodeNumbersR);
    int numOfConn;
    int node1,node2, sign;
    
    map<int,double> nodeVal; // saves the values of the nodes by their nodenumber - needed for sign calculation
    
    for(int i=0; i<numberOfNodes; ++i)
    {
        nodeVal[INTEGER(nodeNumbersR)[i]] = REAL(startValue)[i];
    }
    
    // go through all the nodes
    for(int i=0; i<numberOfNodes; ++i)
    {
        connOneNode = VECTOR_ELT(connR,i); // get the connection of the current node
        numOfConn = LENGTH(connOneNode); // how many are there
        node1 = INTEGER(nodeNumbersR)[i]; // number of node currently working on
        for(int j=0; j<numOfConn; ++j) // go through all connections of this node
        {
            node2 = INTEGER(connOneNode)[j];
            // only check nodes that have a larger number (to avoid doubles and nodes to itself)
            if(node2>node1)
            {
                sign = signum(nodeVal[node1]-nodeVal[node2]);
                if(sign==0)
                {
                    sign = 1; // flow between unconnected nodes should always be -1 or 1
                }
                addEdge(node1,node2,sign);
            }
        }
    }
}

PenaltyGraph::PenaltyGraph(vector<int> nodeNumbers, vector<list<int> > conn, vector<double> startValues)
{
    int numberOfNodes = nodeNumbers.size();
    list<int> connOneNode;
    int numOfConn;
    int node1,node2, sign;
    
    map<int,double> nodeVal; // saves the values of the nodes by their nodenumber - needed for sign calculation
    
    for(int i=0; i<numberOfNodes; ++i)
    {
        nodeVal[nodeNumbers[i]] = startValues[i];
    }
    
    // go through all the nodes
    for(int i=0; i<numberOfNodes; ++i)
    {
        connOneNode = conn[i]; // get the connection of the current node
        numOfConn = connOneNode.size(); // how many are there
        node1 = nodeNumbers[i]; // number of node currently working on
        for(int j=0; j<numOfConn; ++j) // go through all connections of this node
        {
            node2 = connOneNode.front();
            connOneNode.pop_front();
            // only check nodes that have a larger number (to avoid doubles and nodes to itself)
            if(node2>node1)
            {
                sign = signum(nodeVal[node1]-nodeVal[node2]);
                addEdge(node1,node2,sign);
            }
        }
    }
}


PenaltyGraph::~PenaltyGraph()
{
    Node::iterator edgeIt; // iterator over the edges
    Nodes::iterator nodeIt; // iterator over the nodes
    
    for(nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
    {
        for(edgeIt = nodeIt->second.begin(); edgeIt != nodeIt->second.end(); ++edgeIt)
        {
            delete edgeIt->second;
        }
    }
}


int PenaltyGraph::getMaxNodeNum()
{
    Nodes::reverse_iterator iter;
    iter=nodes.rbegin(); // points to the end of the map which is the largest value
    return(iter->first);
}


set<int> PenaltyGraph::allNodes()
{
    set<int> all;
    Nodes::iterator nodeIt;
    
    for(nodeIt=nodes.begin(); nodeIt!=nodes.end(); ++nodeIt)
    {
        all.insert(nodeIt->first);
    }
    return(all);
}





/*
int main(int argc, char** argv)
{
    PenaltyGraph g;
    
    g.addEdge(1,2,1);
    g.addEdge(1,3,1);
    g.addEdge(1,4,1);
    g.addEdge(2,5,1);
    g.addEdge(2,6,1);
    g.addEdge(2,7,1);
    g.addEdge(3,5,1);
    g.addEdge(3,6,1);
    g.addEdge(3,7,1);
    g.addEdge(4,5,1);
    g.addEdge(5,8,1);
    g.addEdge(6,8,1);
    g.addEdge(7,8,1);
    
    g.printGraph(cout);
    
    set<int> subNodes;
    subNodes.insert(5);
    subNodes.insert(6);
    subNodes.insert(7);
    
    MaxFlowGraph *m = g.subGraph(subNodes);
    m->printGraph(cout);
}
    
    vector<MaxFlowGraph*> m(9);
    MaxFlowGraph* pm;
    Groups grps(8);
    
    for(int i=1; i<=8; ++i)
    {
        subNodes.clear();
        subNodes.insert(i);
        m[i-1]=g.subGraph(subNodes);
        grps.addNewGroup(0,0,m[i-1],true);
    }

    grps.printGroups();
    
    // merge two groups
    subNodes.clear();
    subNodes.insert(2);
    subNodes.insert(3);
    pm=g.subGraph(subNodes);
    grps.mergeGroups(1,2,1,pm);
    
    grps.printGroups();
    

    MaxFlowGraph* pm1;
    MaxFlowGraph* pm2;
    subNodes.clear();
    subNodes.insert(2);
    pm1=g.subGraph(subNodes);
    subNodes.clear();
    subNodes.insert(3);
    pm2=g.subGraph(subNodes);
    
    pair<int,int> foo;
    foo=grps.splitGroup(8,2,pm1,pm2);
    cout << foo.first << foo.second << endl;
    
    grps.printGroups();
    
    set<int> bar;
    subNodes.clear();
    subNodes.insert(1);
    subNodes.insert(2);
    subNodes.insert(3);
    bar = grps.nodesToGroups(subNodes);
    set<int>::iterator iter;
    cout << "Translated into Groups: ";
    for(iter=bar.begin(); iter!=bar.end(); ++iter)
    {
        cout << *iter << " ";
    }
    cout << endl;
    
    vector<double> bar2(4);
    vector<double>::iterator iter2;
    bar2[0]=.5;
    bar2[1]=1;
    bar2[2]=1.5;
    bar2[3]=2;
    bar2=grps.nodeSolution(2,bar2);
    for(int i=0; i<4; ++i)
    {
        cout << bar2[i] << " ";
    }
    cout << endl;
    
    for(int i=0; i<1.0e9; ++i)
    {
    }
    
    
    // manipulate the nodes so that they have the capacities we want
    m->nodes[source][2]->capacity=15;
    m->nodes[2][source]->capacity=0;
    m->nodes[source][3]->capacity=10;
    m->nodes[3][source]->capacity=0;
    m->nodes[source][4]->capacity=12; //12
    m->nodes[4][source]->capacity=0;
    m->nodes[2][5]->capacity=5;
    m->nodes[5][2]->capacity=0;
    m->nodes[2][6]->capacity=5;
    m->nodes[6][2]->capacity=0;
    m->nodes[2][7]->capacity=5;
    m->nodes[7][2]->capacity=0;
    m->nodes[3][5]->capacity=6;
    m->nodes[5][3]->capacity=0;
    m->nodes[3][6]->capacity=6;
    m->nodes[6][3]->capacity=0;
    m->nodes[3][7]->capacity=6;
    m->nodes[7][3]->capacity=0;
    m->nodes[4][5]->capacity=12;
    m->nodes[5][4]->capacity=0;
    m->nodes[5][sink]->capacity=10;
    m->nodes[sink][5]->capacity=0;
    m->nodes[6][sink]->capacity=15;
    m->nodes[sink][6]->capacity=0;
    m->nodes[7][sink]->capacity=15;
    m->nodes[sink][7]->capacity=0;

    
//    m->setFlowTo0();
}
*/
