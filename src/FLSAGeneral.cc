#include <fstream>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "FLSAGeneral.h"
#include "MaxFlowGraph.h"

void FLSAGeneral::initializeGroups(SEXP connList, SEXP startValues)
{
    // set up the groups
    SEXP nodeNumbersR = VECTOR_ELT(connList,0);
    int numOfNodes = LENGTH(nodeNumbersR);
    MaxFlowGraph* m;
    set<int> subNodes;
    
    // print progress report if necessary
    if(showProgress)
    {
        Rprintf("Started initializing the Groups\n");
    }
    
    for(int i=0; i<numOfNodes; ++i)
    {
        // the node in the subgraph
        subNodes.clear();
        subNodes.insert(INTEGER(nodeNumbersR)[i]);
        // the subgraph
        m=graph.subGraph(subNodes);
        
        // add it as a new group at lambda=0
        groups.addNewGroup(0,REAL(startValues)[i], m, true);
    }
    if(showProgress)
    {
        Rprintf("Finished initializing the Groups\n");
    }
}




void FLSAGeneral::initializeScheduler()
{
    set<int> allNodes = graph.allNodes(); // get all nodes that are in the penaltygraph
    set<int> curNode, connNodes,bar, connGroups;
    groupItem curGroupItem, connGroupItem;
    int curGroupNum;
    
    set<int>::iterator setIter, setIter2;
    // go through all nodes
    
    if(showProgress)
    {
        Rprintf("Started initializing the Scheduler\n");
    }

    for(setIter=allNodes.begin(); setIter!=allNodes.end(); ++setIter)
    {
        curNode.clear();
        curNode.insert(*setIter);
        bar = groups.nodesToGroups(curNode);
        curGroupNum = *(bar.begin()); // group the current node is in
        
        // which groups is the current node connected to
        connGroups = groups.nodesToGroups(graph.connectedTo(curNode));
        // erase all with number less than or equal to the current group to avoid doubles
        connGroups.erase(connGroups.begin(), connGroups.upper_bound(curGroupNum));
        
        // since all groups consist of only one member, tension updating is not needed;
        
        // schedule all the events of the current group and connected ones
        scheduleMergeEvents(curGroupNum, connGroups);
    }

    if(showProgress)
    {
        Rprintf("Finished initializing the Scheduler\n");
    }

}



double FLSAGeneral::calcHitTime(groupItem grp1, groupItem grp2)
{
    double lambdaMax = Max(grp1.lambda, grp2.lambda);
    double lhs, slopeRhs, offset;
    lhs = grp1.mu-grp2.mu+grp1.deriv*(lambdaMax-grp1.lambda)-grp2.deriv*(lambdaMax-grp2.lambda);
    slopeRhs = grp2.deriv-grp1.deriv;
    // rhs = (grp2.deriv-grp1.deriv)*(lambda-lambdaMax)
    if(showProgress)
    {
        Rprintf("LHS: %f RHS: %f\n",lhs,slopeRhs);
        Rprintf("Group 1: Lambda: %f Deriv: %f Size: %d\n", grp1.lambda, grp1.deriv, grp1.size);
        Rprintf("Group 2: Lambda: %f Deriv: %f Size: %d\n", grp2.lambda, grp2.deriv, grp2.size);
    }

    if(RelDif(lhs,0)<tolerance) // lhs is 0, so hitting time is lambdaMax as groups on top of each other
    {
        // check that the groups will actually hit and are not already past each other
        // for this it is necessary that the two groups are connected (otherwise throws an exception)
        int flowSign = graph.flowSignBetweenGroups(grp1.m->allNodes(), grp2.m->allNodes());
        int rhsSign = signum(slopeRhs);
        if(showProgress)
        {
            Rprintf("FlowSign: %d; rhsSign: %d", flowSign, rhsSign);
        }

        if((flowSign==0) || (rhsSign==0)) // groups on top of each other and don't move; merge them
        {
            return(lambdaMax);
        }
        else if(flowSign == rhsSign) // will meet immediately
        {
            return(lambdaMax);
        }
        else // flowSign != rhsSign
        {
            return(neverHit);
        }
    }
    else // lhs is not 0
    {
        if(RelDif(slopeRhs,0)<tolerance) // slope is 0
        {
            return(neverHit);
        }
        else
        {
            offset = lhs/slopeRhs;
            if(offset < -tolerance) // only hits in the past
            {
                return(neverHit);
            }
            else
            {
                return(offset + lambdaMax);
            }
        }
    }
}


// initialize the group class to have sufficient space such that the largest nodeNumber can still be stored
FLSAGeneral::FLSAGeneral(int highestNodeNum, SEXP connList, SEXP startValues, SEXP splitCheckSize, SEXP verbose, SEXP tol, SEXP maxGrpNum, double highestLambda):groups(highestNodeNum+1),graph(connList, startValues)
{
    maxLambda=highestLambda; // maximum lambda up to which to compute
    maxSizeForSplitCheck = INTEGER(splitCheckSize)[0]; // above this value, groups won't be checked for splits
    showProgress = LOGICAL(verbose)[0];
    tolerance = REAL(tol)[0];
    maxGroupNumber = INTEGER(maxGrpNum)[0];
    // set up the initial groups
    initializeGroups(connList, startValues);
    
    // insert the events of the first groups into the scheduler
    initializeScheduler();
    
    runAlgorithm();
    
}

void FLSAGeneral::runAlgorithm()
{
    pair<double, scheduleEvent> schedEv;
    
    // run as long as there is still an event in the scheduler left
    while(!scheduler.empty() && groups.size() < maxGroupNumber)
    {
        R_CheckUserInterrupt();
        schedEv = scheduler.getNextEvent();
        
        // if lambda in here is greater than lambdaMax, finish executing the algorithm
        if(schedEv.first>maxLambda)
        {
            return;
        }
        
        // possible event types are 'M' for merge 'T' for tension reevaluation (might lead to splitting)
        if(schedEv.second.type=='M') // do the merging
        {
            doMerging(schedEv.first, schedEv.second.grp1, schedEv.second.grp2);
        }
        else if(schedEv.second.type=='T') // reevaluate tension changes
        {
            doTension(schedEv.first, schedEv.second.grp1, true);
        }
        else
        {
            throw("wrong type in schedule event");
        }
    }
    if(groups.size() >= maxGroupNumber) 
    {
        error("Number of groups too large. Try increasing the tolerance!\n");
    }
}


void FLSAGeneral::scheduleMergeEvents(int grpNum, const set<int>& connGroups)
{
    set<int>::const_iterator setIter;
    groupItem curGroupItem = groups.getGroup(grpNum);
    // schedule new events for these groups
    for(setIter=connGroups.begin(); setIter!=connGroups.end(); ++setIter)
    {   // get the group it is connected to, calculate the hit time and schedule
        groupItem connGroupItem = groups.getGroup(*setIter);
        double hitTime = calcHitTime(curGroupItem, connGroupItem);
        if(hitTime!=neverHit) // a hit does occur; schedule it
        {
            scheduleEvent schedEv;
            schedEv.type='M';
            schedEv.grp1=grpNum;
            schedEv.grp2= *setIter;
            scheduler.insertEvent(hitTime, schedEv);
        }
    }
}



void FLSAGeneral::doMerging(double lambda, int grp1, int grp2)
{
//    outStream << "Just started in doMerging" << endl;
    // check if both groups are still active
    if(groups.isActive(grp1) && groups.isActive(grp2))
    {
        set<int> grp1Nodes, grp2Nodes, connGroups;
        groupItem foo;
        MaxFlowGraph* m;
        int newGroupNum;
        set<int>::iterator setIter;

        // get the nodes in the two groups
        foo = groups.getGroup(grp1);
        grp1Nodes = foo.m->allNodes();
        foo = groups.getGroup(grp2);
        grp2Nodes = foo.m->allNodes();
        
        // print a status message
        if(showProgress)
        {
            Rprintf("Lambda: %f Action: M Groups: %d, %d Sizes: %d, %d\n",lambda, grp1, grp2, grp1Nodes.size(), grp2Nodes.size());
        }
        
        // merge the nodes and get the new subgraph and set up the new group
        grp1Nodes.insert(grp2Nodes.begin(), grp2Nodes.end());
        m = graph.subGraph(grp1Nodes);
        newGroupNum = groups.mergeGroups(grp1, grp2, lambda,m);
        // which groups is that new graph connected to
        connGroups = groups.nodesToGroups(graph.connectedTo(grp1Nodes));
        // schedule the new events
        scheduleMergeEvents(newGroupNum, connGroups);
        
        // merging also makes an tension update necessary
        doTension(lambda, newGroupNum, false);
    
    }
}




void FLSAGeneral::doTension(double lambda, int grp, bool update)
{
    double hitTime;
    // check that the group is still active
    if(groups.isActive(grp))
    {
         
        groupItem curGroupItem = groups.getGroup(grp);
        // do not do any tension updates if the group size is too large
        if(showProgress)
        {
            Rprintf("Lambda: %f Action: T Group: %d Size: %d\n",lambda, grp, curGroupItem.m->size());
        }
       if(curGroupItem.m->size()<=maxSizeForSplitCheck)
        {
            // for this group, first initiate a tension change calculation
            if(update)
            {
                hitTime = curGroupItem.m->calcTensionChangeUpdate(lambda);
            }
            else
            {
                hitTime = curGroupItem.m->calcTensionChangeProportional(lambda);
            }
 
            // check the different possibilities for the hitTime
            if(hitTime==neverSplit) // do not schedule anything
            {
                return;
            }
            else if(hitTime==splitNow) // make the split and schedule the new events
            {
                split(lambda, grp);
                return;
            }
            else // schedule the new tension event
            {
                scheduleEvent schedEvNew;
                schedEvNew.type='T';
                schedEvNew.grp1=grp;
                scheduler.insertEvent(hitTime, schedEvNew);
                return;
            }
        }
    } 
}



void FLSAGeneral::split(double lambda, int grp)
{
//    outStream << "Just started in Split" << endl;
    // get the group that is to be split
    if(showProgress)
    {
        Rprintf("Lambda: %f Action: Split Group: %d\n",lambda, grp);
    }

    
    groupItem groupToSplit = groups.getGroup(grp);
    
    // save the nodes of the first splitted graph and the second splitted graph
    set<int> splitNodes1, splitNodes2;
    splitNodes1 = groupToSplit.m->reachableFromSource();
    splitNodes2 = groupToSplit.m->getComplement(splitNodes1);
    
    // generate the two new subgraphs
    MaxFlowGraph *mSplit1, *mSplit2;
    mSplit1 = graph.subGraph(splitNodes1);
    mSplit2 = graph.subGraph(splitNodes2);
    // do the splitting
    pair<int, int> newGroupNums = groups.splitGroup(grp,lambda,mSplit1, mSplit2);
    
    // for each of the new groups find the groups they are connected to
    set<int> connGroups1 = groups.nodesToGroups(graph.connectedTo(splitNodes1));
    set<int> connGroups2 = groups.nodesToGroups(graph.connectedTo(splitNodes2));
    
    // the splitted groups should not be merged again immediately; so delete them from the
    // list of groups they are connected to
    connGroups1.erase(newGroupNums.second);
    connGroups2.erase(newGroupNums.first);
    
    // Schedule all possible future merge events
    scheduleMergeEvents(newGroupNums.first, connGroups1);
    scheduleMergeEvents(newGroupNums.second, connGroups2);
    
    // update the tensions in the groups
    doTension(lambda, newGroupNums.first, false);
    doTension(lambda, newGroupNums.second, false);
}



SEXP FLSAGeneral::solution(SEXP nodes, SEXP lambdas)
{
    return(groups.solution(nodes,lambdas));
}


SEXP FLSAGeneral::solutionGraph()
{
    return(groups.getSolutionObject());
}



/* FLSAGeneral::FLSAGeneral(int highestNodeNum, vector<int> nodeNumbers, vector<list<int> > conn, vector<double> startValues, int splitCheckSize):groups(highestNodeNum+1),graph(nodeNumbers, conn, startValues),outStream("TestOutput.txt", ios::out)
{
    maxSizeForSplitCheck = splitCheckSize;
//    cout << "Before the opening of the file" << endl;
    // set up the initial groups
    initializeGroups(nodeNumbers,  startValues);
//    groups.printGroups(outStream);
    // insert the events of the first groups into the scheduler
    initializeScheduler();
    
//    scheduler.printSchedule(outStream);
    
    // run the algorithm
    runAlgorithm();
    
//    groups.printGroups(outStream);
}

void FLSAGeneral::initializeGroups(vector<int> nodeNumbers, vector<double> startValues)
{
    // set up the groups
    int numOfNodes = nodeNumbers.size();
    MaxFlowGraph* m;
    set<int> subNodes;
    
    for(int i=0; i<numOfNodes; ++i)
    {
        // the node in the subgraph
        subNodes.clear();
        subNodes.insert(nodeNumbers[i]);
        // the subgraph
        m=graph.subGraph(subNodes);
        
        // add it as a new group at lambda=0
        groups.addNewGroup(0,startValues[i], m, true);
    }
}

*/
