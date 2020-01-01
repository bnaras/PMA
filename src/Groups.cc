#include <iostream>
#include "Groups.h"
#include "GeneralFunctions.h"

void Groups::inactivateGroup(int grp, double lambda)
{
  if(((unsigned int) grp) >= groups.size()) // nothing to inactivate as this group does not exist
    {
        return;
    }
    
    if(!groups[grp].active) // nothing to inactivate as group is already inactive
    {
        return;
    }
    
    groups[grp].active=false;
    groups[grp].endLambda=lambda;
    delete(groups[grp].m);
}


void Groups::updateNodeMap(set<int> nodes, int grp, bool initial)
{
    set<int>::iterator iter;
    
    
    for(iter=nodes.begin(); iter!=nodes.end(); ++iter)
    {
        nodeMap[*iter]=grp;
    }
    if(initial)
    {
        for(iter=nodes.begin(); iter!=nodes.end(); ++iter)
        {
            initialNodeMap[*iter]=grp;
        }
    }
}

SEXP Groups::solutionObjectInit()
{
    SEXP solObj;
    PROTECT(solObj = allocVector(VECSXP,11));
    // give the 8 elements names 
    SEXP names;
    PROTECT(names = allocVector(STRSXP,11));
    SET_STRING_ELT(names,0,mkChar("Number"));
    SET_STRING_ELT(names,1,mkChar("BeginLambda"));
    SET_STRING_ELT(names,2,mkChar("EndLambda"));
    SET_STRING_ELT(names,3,mkChar("Mu"));
    SET_STRING_ELT(names,4,mkChar("Derivative"));
    SET_STRING_ELT(names,5,mkChar("Action"));
    SET_STRING_ELT(names,6,mkChar("Group1"));
    SET_STRING_ELT(names,7,mkChar("Group2"));
    SET_STRING_ELT(names,8,mkChar("Group1Nodes"));
    SET_STRING_ELT(names,9,mkChar("Size"));
    SET_STRING_ELT(names,10,mkChar("InitialNodeMap"));
    setAttrib(solObj, R_NamesSymbol, names);
    UNPROTECT(1);
    
    // give the object a class
    SEXP className;
    PROTECT(className=allocVector(STRSXP,1));
    SET_STRING_ELT(className,0,mkChar("FLSAGeneral"));
    classgets(solObj,className);
    UNPROTECT(1);
    
    // install the right vectors in the list
    int numGroups = groups.size();
    SET_VECTOR_ELT(solObj,0,allocVector(INTSXP,numGroups));
    SET_VECTOR_ELT(solObj,1,allocVector(REALSXP,numGroups));
    SET_VECTOR_ELT(solObj,2,allocVector(REALSXP,numGroups));
    SET_VECTOR_ELT(solObj,3,allocVector(REALSXP,numGroups));
    SET_VECTOR_ELT(solObj,4,allocVector(REALSXP,numGroups));
    SET_VECTOR_ELT(solObj,5,allocVector(INTSXP,numGroups));
    SET_VECTOR_ELT(solObj,6,allocVector(INTSXP,numGroups));
    SET_VECTOR_ELT(solObj,7,allocVector(INTSXP,numGroups));
    SET_VECTOR_ELT(solObj,8,allocVector(VECSXP,numGroups));
    SET_VECTOR_ELT(solObj,9,allocVector(INTSXP,numGroups));
    SET_VECTOR_ELT(solObj,10,allocVector(INTSXP,initialNodeMap.size()));
    UNPROTECT(1);
    return(solObj);
}

SEXP Groups::getSolutionObject()
{
    // allocate a list with 2 elements
    SEXP solObj;
    PROTECT(solObj = solutionObjectInit());
    
    // now go through all groups and save the relevant data
    groupItem grpObj;
    int numGroups = groups.size();
    for(int grp =0; grp<numGroups; ++grp)
    {
        grpObj = getGroup(grp);
        INTEGER(VECTOR_ELT(solObj,0))[grp]=grp;
        REAL(VECTOR_ELT(solObj,1))[grp]=grpObj.lambda;
        REAL(VECTOR_ELT(solObj,2))[grp]=grpObj.endLambda;
        REAL(VECTOR_ELT(solObj,3))[grp]=grpObj.mu;
        REAL(VECTOR_ELT(solObj,4))[grp]=grpObj.deriv;
        if(grpObj.action=='M')
        {
            INTEGER(VECTOR_ELT(solObj,5))[grp]=0;
        }
        else // a split event
        {
            INTEGER(VECTOR_ELT(solObj,5))[grp]=1;
            // also set the 8th vector element
            int numSplitNodes = grpObj.splitNodes.size();
            set<int>::iterator setIt;
            int i;
            SET_VECTOR_ELT(VECTOR_ELT(solObj,8),grp,allocVector(INTSXP,numSplitNodes));
            for(i=0, setIt = grpObj.splitNodes.begin(); i<numSplitNodes;++setIt, ++i)
            {
                INTEGER(VECTOR_ELT(VECTOR_ELT(solObj,8),grp))[i]=*setIt;
            }
        }
        INTEGER(VECTOR_ELT(solObj,6))[grp]=grpObj.grp1;
        INTEGER(VECTOR_ELT(solObj,7))[grp]=grpObj.grp2;
        INTEGER(VECTOR_ELT(solObj,9))[grp]=grpObj.size;
    }
    // transfer the initial mapping of the nodes
    for(unsigned int i=0; i<initialNodeMap.size(); ++i)
    {
        INTEGER(VECTOR_ELT(solObj,10))[i]=initialNodeMap[i];
    }
    
    UNPROTECT(1);
    return(solObj);
}

Groups::Groups(int numNodes)
{
    nodeMap.assign(numNodes, noGroup);
    initialNodeMap.assign(numNodes, noGroup);
}

Groups::Groups(SEXP solution)
{
    int numGroups = LENGTH(VECTOR_ELT(solution,0));
    groupItem nullGroup = {0};
    
    groups.assign(numGroups, nullGroup);
    
    for(int i=0; i<numGroups; ++i)
    {
        // write all the basic data into the right variables
        groups[i].lambda = REAL(VECTOR_ELT(solution,1))[i];
        groups[i].endLambda = REAL(VECTOR_ELT(solution,2))[i];
        groups[i].mu = REAL(VECTOR_ELT(solution,3))[i];
        groups[i].deriv = REAL(VECTOR_ELT(solution,4))[i];
        if(INTEGER(VECTOR_ELT(solution,5))[i]==1)
        {
            groups[i].action = 'S';
        }
        else
        {
            groups[i].action = 'M';
        }
        groups[i].grp1 = INTEGER(VECTOR_ELT(solution,6))[i];
        groups[i].grp2 = INTEGER(VECTOR_ELT(solution,7))[i];
        groups[i].size = INTEGER(VECTOR_ELT(solution,9))[i];
        
        // write the splitNodes back into the data
        if(groups[i].action=='S')
        {
            int numSplitNodes = LENGTH(VECTOR_ELT(VECTOR_ELT(solution,8),i));
            for(int j=0; j<numSplitNodes; ++j)
            {
                groups[i].splitNodes.insert(INTEGER(VECTOR_ELT(VECTOR_ELT(solution,8),i))[j]);
            }
        }
    }
        // write the inital nodeMap
    int initialNodeMapSize = LENGTH(VECTOR_ELT(solution,10));
    initialNodeMap.assign(initialNodeMapSize,noGroup);
    for(int j=0; j<initialNodeMapSize; ++j)
    {
        initialNodeMap[j] = INTEGER(VECTOR_ELT(solution,10))[j];
    }

}


int Groups::addNewGroup(double lambda, double mu, MaxFlowGraph* m, bool initial)
{
    groupItem g;
    // fill out the group item
    g.lambda = lambda;
    g.mu=mu;
    g.active=true;
    g.deriv = m->getGroupDeriv();
    g.m = m;
    g.endLambda=infinite;
    g.size = m->size();
    
    // include g in the vector
    int groupNum = groups.size();
    groups.push_back(g);
    
    // update the node - group relationship
    updateNodeMap(m->allNodes(), groupNum, initial);
    return(groupNum);
}




int Groups::mergeGroups(int grp1, int grp2, double lambda, MaxFlowGraph* m)
{
    // inactivate the 2 old groups
    inactivateGroup(grp1,lambda);
    inactivateGroup(grp2,lambda);
    
    // calculate the new mu
    double mu = groups[grp1].mu + groups[grp1].deriv*(lambda-groups[grp1].lambda);
    
    // make the new group
    int newGrp = addNewGroup(lambda, mu, m);
    
    // set the action in the old groups correctly
    groups[grp1].action='M';
    groups[grp1].grp1=newGrp;
    groups[grp1].grp2=0;
    groups[grp2].action='M';
    groups[grp2].grp1=newGrp;
    groups[grp2].grp2=0;
    
    // return the number of the new group
    return(newGrp);
}




pair<int, int> Groups::splitGroup(int grp, double lambda, MaxFlowGraph* m1, MaxFlowGraph* m2)
{
    // inactivate the old group
    inactivateGroup(grp,lambda);
    
    // calculate the new mu
    double mu = groups[grp].mu + groups[grp].deriv*(lambda-groups[grp].lambda);
    
    // make the 2 new groups
    pair<int,int> newGrps;
    newGrps.first = addNewGroup(lambda, mu, m1);
    newGrps.second = addNewGroup(lambda, mu, m2);
    
    // set the action in the old correctly
    groups[grp].action = 'S';
    groups[grp].grp1=newGrps.first;
    groups[grp].grp2=newGrps.second;
    groups[grp].splitNodes = m1->allNodes();
    
    // return the numbers of the two new groups
    return(newGrps);
}




set<int> Groups::nodesToGroups(const set<int>& nodes)
{
    set<int> grpSet;
    set<int>::const_iterator iter;
    
    // go through all the nodes in the list and translate them into the currently active group
    for(iter=nodes.begin(); iter!=nodes.end(); ++iter)
    {
        if(nodeMap[*iter]!=noGroup) 
        {
            grpSet.insert(nodeMap[*iter]);
        }
        else // every node has to belong to a group
        {
            throw("Asked for node that does not belong to a group in 'nodesToGroup'");
        }
    }
    
    // return the set of groups
    return(grpSet);
}






vector<double> Groups::nodeSolution(int node,const vector<double>& lambdas)
{
    // vector to save the results
    vector<double> solution(lambdas.size(), 0);
    
    // initialize index and group counters
    int numLambdas = lambdas.size();
    int curLambdaIndex = 0;
    int curGrp=initialNodeMap[node];
    
    //check that the node was ever assigned (checks for possible input errors with respect to the nodes)
    if(curGrp==noGroup)
    {
        throw("Node asked for in Groups::nodeSolution was never assigned a group.");
    }
    
    while(curLambdaIndex< numLambdas) // keep going as long as more lambdas are left
    {
        // check if in the group that covers the current lambda
        if(groups[curGrp].endLambda>=lambdas[curLambdaIndex])
        {
            solution[curLambdaIndex] = groups[curGrp].mu + groups[curGrp].deriv*(lambdas[curLambdaIndex]-groups[curGrp].lambda);
            ++curLambdaIndex;
        }
        else // go up to the next group
        {
            // if it is a merge, just go up to the next group
            if(groups[curGrp].action=='M')
            {
                curGrp=groups[curGrp].grp1;
            }
            else if(groups[curGrp].action=='S') // if it is a split, decide the correct next group
            {
                // check if node goes into the first or second new grp
                if(groups[curGrp].splitNodes.count(node)>0) // first group
                {
                    curGrp = groups[curGrp].grp1;
                }
                else // second group
                {
                    curGrp = groups[curGrp].grp2;
                }
            }
            else
            {
                throw("Unspecified action type in nodeSolution");
            }
        }
    }
    return(solution);
}



SEXP Groups::solution(SEXP nodes, SEXP lambdas)
{
    int nodesLen = LENGTH(nodes);
    int lambdasLen = LENGTH(lambdas);
    // transfer the lambdas into a vector
    vector<double> lambdasVec(lambdasLen);
    for(int i=0; i<lambdasLen; ++i)
    {
        lambdasVec[i]=REAL(lambdas)[i];
    }
    // get a vector of the right size in R
    // the nodes will be in the columns, the lambdas in the rows
    SEXP sol;
    PROTECT(sol=allocMatrix(REALSXP,lambdasLen, nodesLen));
    double *solMat = REAL(sol);
    // go through all the nodes
    int counter=0;
    for(int i =0; i<nodesLen; ++i)
    {
        vector<double> nodeSol = nodeSolution(INTEGER(nodes)[i],lambdasVec);
        // copy the solution into the R vector
        for(int j=0; j<lambdasLen; ++j)
        {
            solMat[counter]=nodeSol[j];
            ++counter;
        }
    }
    
    // the names for the matrix will be the numbers of the nodes and lamba
    SEXP dimnames;
    PROTECT(dimnames = allocVector(VECSXP,2));
    SET_VECTOR_ELT(dimnames,0,lambdas);
    SET_VECTOR_ELT(dimnames,1,nodes);
    setAttrib(sol,R_DimNamesSymbol, dimnames);
    UNPROTECT(2);
    return(sol);
}


void Groups::printGroups(ostream& outStream)
{
    //print the initial mapping between the nodes and the groups
    outStream << "Initial mapping of the nodes:" << endl;
    for(unsigned int i=0; i<initialNodeMap.size(); ++i)
    {
        outStream << "Node: " << i << " Group: " << initialNodeMap[i] << endl;
    }
    // print the current mapping between the nodes and the groups
    outStream << "Current mapping of the nodes:" << endl;
    for(unsigned int i=0; i<nodeMap.size(); ++i)
    {
        outStream << "Node: " << i << " Group: " << nodeMap[i] << endl;
    }
    

    for(unsigned int i=0; i<groups.size(); ++i)
    {
        outStream << "-------------------------------------------------------" << endl;
        outStream << "Group Number: " << i << endl;
        outStream << "Lambda: " << groups[i].lambda << " Mu: " << groups[i].mu << " Deriv: " << groups[i].deriv << " EndLambda: " << groups[i].endLambda << endl;
        outStream << "Active: " << groups[i].active << " Action: " << groups[i].action << endl;
        outStream << "Group 1: " << groups[i].grp1 << " Group 2: " << groups[i].grp2 << endl;
        outStream << "Split Nodes: ";
        // print all nodes in groups[i].splitNodes
        set<int>::iterator iter;
        for(iter = groups[i].splitNodes.begin(); iter!=groups[i].splitNodes.end(); ++iter)
        {
            outStream << *iter << " ";
        }
        outStream << endl << "MaxFlowGraph" << endl;
        if(groups[i].active)
        {
            //groups[i].m->printGraph(outStream);
        }
        
        outStream << "---------------------------------------------------------------------" << endl;
    }
}
