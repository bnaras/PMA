#ifndef _GROUPS_
#define _GROUPS_


#include <set>
#include <vector>
#include "MaxFlowGraph.h"
#include <R.h>
#include <Rinternals.h>


using namespace std;

const int noGroup = -1;
// define the basic data available for each group
struct groupItem {
    double lambda; // when the group was created
    double mu; // value for the whole group
    double deriv; // derivative of mu
    double endLambda; // value of lambda when the groups stops to exist
    bool active; // is the node still active?
    char action; // which action to take at the end of the group life 'M' for merge 'S' for split
    int grp1; // needed for splitting and merging
    int grp2; // only needed when splitting
    set<int> splitNodes; // used when splitting; contains the nodes that are split into group1; rest into group2
    int size; // stores the size of the group
    
    MaxFlowGraph* m; // maxflowgraph that contains the tension information
};



class Groups
{
    vector<groupItem> groups; // vector in which all the groups will be saved
    vector<int> nodeMap; // stores the mapping of the nodes to the currently active groups;
    vector<int> initialNodeMap; // stores the mapping of nodes to the groups they belong to at lambda=0;
    
    // helper function that inactivates groups (deletes the maxflowgraph and other elements that are not needed)
    void inactivateGroup(int grp, double lambda);
    // update the group to which the nodes belong
    void updateNodeMap(set<int> nodes, int grp, bool initial=false);
    
    // makes the setup for the return of the solution object
    SEXP solutionObjectInit();
    // function that, given a nodenumber and a vector of lambdas, returns a vector of equal 
    // length with the solution
    vector<double> nodeSolution(int node, const vector<double>& lambdas);
    
public:
    // a constructor that is based on the R object that is returned as the results
    Groups(SEXP solution);

    // a function that returns the complete solution as an R object
    SEXP getSolutionObject();

    // constructor needs the number of nodes in the problem
    Groups(int numNodes);

    // add a new group to the data that did not exist before (can only consist of one node)
    // return value is the number of the group (can not be user specified, as this is
    // what the class mainly controls)
    int addNewGroup(double lambda, double mu, MaxFlowGraph* m, bool initial=false);
    
    // merge two given groups into a new one; return value is the number of the new group
    // maxflowgraphs in old groups will be deleted as the group is set inactive
    int mergeGroups(int grp1, int grp2, double lambda, MaxFlowGraph* m);
    
    // split up a group into 2 new groups; return value is a pair of the new group numbers
    // the first of the pair corresponds to m1, the second to m2
    pair<int, int> splitGroup(int grp, double lambda, MaxFlowGraph* m1, MaxFlowGraph* m2);
    
    // return a group; it is not being checked that this group actually exists
    inline groupItem getGroup(int grp){ return(groups[grp]);};
    
    // return if a group is still active
    bool isActive(int grp) {return(groups[grp].active);};
    
    // given a set of nodes, return a set of grp they currently belong to
    // this can be used to find out to which groups a group is connected to
    set<int> nodesToGroups(const set<int>& nodes);
    
    // given a set of nodes and a list of sorted lambdas (increasing), the function
    // returns a matrix 
    SEXP solution(SEXP nodes, SEXP lambdas);

    // print all the data in all the groups; will be useful for debugging
    void printGroups(ostream& outStream);

    inline int size() {return(groups.size());};
};


#endif
