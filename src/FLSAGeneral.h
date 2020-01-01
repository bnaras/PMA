/*********************************************************************
***
*** This class implements the generalized Fused Lasso Signal Approximator Algorithm
*** The class itself will be given teh necessary R objects in its constructor
*** Then the whole solution path will be calculated
*** The only used accessible function will return the solution tree or the solution values
*** itself on request as R objects
***
*********************************************************************/
#ifndef _FLSAGENERAL_
#define _FLSAGENERAL_

#include "PenaltyGraph.h"
#include "MaxFlowGraph.h"
#include "Scheduler.h"
#include "Groups.h"
#include <set>


const double neverHit=-1;

class FLSAGeneral
{
    Groups groups;
    PenaltyGraph graph;
    Scheduler scheduler;
    int maxSizeForSplitCheck;
    bool showProgress;
    double maxLambda;
    double tolerance;
    double maxGroupNumber;
    
    // helper function for the constructor that performs the initialization of the groups
    void initializeGroups(SEXP connList, SEXP startValues);

    // helper function for the constructor that performs the initialization of the scheduler
    void initializeScheduler();
    
    // given two group items, calculates if and when the 2 groups will hit (only in the future)
    double calcHitTime(groupItem grp1, groupItem grp2);
    
    // run the algorithm until there is only one group left
    void runAlgorithm();
    
    // given a number of the new group and a list of groups it is connected to, schedule the new events
    void scheduleMergeEvents(int grpNum,const set<int>& connGroups);
    
    // given a scheduled event, do the merging; does not check anymore that it actually is a merging event
    void doMerging(double lambda, int grp1, int grp2);
    
    // given a scheduled event, do the tension updating; 
    // the return value is the return value of tensionChange in the maxflowgraph class and has to be processed;
    // standard processing for the return value with splitting or scheduling a tension event or doing nothing will be
    // done by the function tensionPostProcessing
    void doTension(double lambda, int grp, bool update); // used for tension after merging or splitting
    // with update = false, for update events with update = true
    
    // function that given lambda and the group number performs the splitting of the group
    // should only be called from function doTension if splitting was requested;
    void split(double lambda, int grp);
    
    
public:

    // constructor 
    FLSAGeneral(int highestNodeNum, SEXP connList, SEXP startValues, SEXP splitCheckSize, SEXP verbose, SEXP tol, SEXP maxGrpNum, double highestLambda=infinite);
    
    // given a set of nodes and a list of sorted lambdas (increasing), the function
    // returns a matrix 
    SEXP solution(SEXP nodes, SEXP lambdas);

    // returns an R object that describes the whole solution and can later be used to extract
    // values for certain lambdas
    SEXP solutionGraph();
};







#endif

/*
FLSAGeneral(int highestNodeNum, vector<int> nodeNumbers, vector<list<int> > conn, vector<double> startValues, int splitCheckSize);
void initializeGroups(vector<int> nodeNumbers, vector<double> startValues);
*/
