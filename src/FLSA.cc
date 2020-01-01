#include <R.h>
#include <Rinternals.h>
#include <vector>
#include <map>
#include "FLSA.h"
#include "GeneralFunctions.h"

using namespace std;




void FLSAClass::addConnection(int grpOne, int grpTwo, double lambda)
{
    double  y1,y2, deriv1, deriv2;
    Connection conn;
    
    y1 = getCurrMu(groupVec[grpOne],lambda);
    y2 = getCurrMu(groupVec[grpTwo],lambda);
    deriv1 =groupVec[grpOne].deriv;
    deriv2 = groupVec[grpTwo].deriv;
   
    conn.first = grpOne;
    conn.second= grpTwo;

    // see if the difference is 0
    if(RelDif(y1,y2) < tolerance)
    {
        // insert into the multimap
        groupMove.insert(make_pair(lambda, conn));
    }
    else // calculate when the two groups would hit
    {
        if(RelDif(deriv1,deriv2)>=tolerance) // will only hit if the derivatives are not equal
        {
            double lambdaCorr = -(y1-y2)/(deriv1-deriv2);
            if(lambdaCorr > 0) // accuracy not necessary here as taken care of before
            groupMove.insert(make_pair(lambdaCorr + lambda, conn));
        }
    }
}


vector<int> FLSAClass::getNeighbours(int grpNum, int exclGrp)
{
    typedef vector<int>::iterator VI;
    vector<int> neighbour, res;
    
    neighbour = groupVec[grpNum].neighbour;
    
    for(VI i=neighbour.begin(); i!=neighbour.end(); ++i)
    {
        if(*i != exclGrp)
        {
            res.push_back(*i);
        }
    }
    return(res);
}

void FLSAClass::updateNeighbours(vector<int> updateGrp, int oldGrp, int newGrp)
{
    for(unsigned int i=0; i< updateGrp.size(); ++i)
    {
        // number of neighbours in the group to update
        int numNeighbours = groupVec[updateGrp[i]].neighbour.size();

        // look at every neighbour and update the numbers from old to new if necessary
        for(int j=0; j< numNeighbours; ++j)
        {
            if(groupVec[updateGrp[i]].neighbour[j]==oldGrp)
            {
                groupVec[updateGrp[i]].neighbour[j] = newGrp;
            }
        }
    }
}


void FLSAClass::deactivateGroup(int grpNum, int newGrp, double lambda)
{
    groupVec[grpNum].active=false;
    groupVec[grpNum].mergeTo=newGrp;
    groupVec[grpNum].mergeLambda=lambda;
}




pair<double, Connection> FLSAClass::getNextConnection()
{
    pair<double, Connection> res;
    typedef map<double, Connection>::iterator MI;
    MI p; // iterator pointing at the first connection
    
    // grab the connection with the lowest lambda value until both
    // groups of the connection are still active
    for(;;) // break when the connection has been found
    {
        p=groupMove.begin();
        // check if any elements are left
        if(p==groupMove.end())
        {
            // no left
            res.first = -1; // signals that nothing could be retrieved (lambda >=0 otherwise)
            break;
        }
        res = *p; // copy it
        groupMove.erase(p); // erase it
        // now check that both groups are still active
        if(groupVec[res.second.first].active && groupVec[res.second.second].active)
        {
            break;
        }
    }
    return(res);
}


FLSAClass::FLSAClass(SEXP yR)
{
    // check if y is ok
    checkInput(yR);
    // calculate parameter length
    int n = LENGTH(yR);
    double* y = REAL(yR);
    
    // save the number of variables
    numVariables = n;

    // there will be a total of 2n-1 groups
    groupVec.resize(2*n-1);
    maxgroup=n-1;
    // initialize the first n elements
    for(int i =0; i< n; i++)
    {
        groupVec[i].active=true;
        groupVec[i].mu=y[i];
        groupVec[i].lambda=0;
        groupVec[i].grpSize = 1;
        groupVec[i].mergeLambda=-1;
        groupVec[i].mergeTo=-1;
        
        if(i==0)
        {
            groupVec[i].neighbour.resize(1);
            groupVec[i].neighbour[0]=1;
            groupVec[i].deriv = signum(y[1]-y[0]);
        }
        else if(i==(n-1))
        {
            groupVec[i].neighbour.resize(1);
            groupVec[i].neighbour[0]=n-2;
            groupVec[i].deriv = signum(y[n-2]-y[n-1]);
        }
        else
        {
            groupVec[i].neighbour.resize(2);
            groupVec[i].neighbour[0]=i-1;
            groupVec[i].neighbour[1]=i+1;
            groupVec[i].deriv = signum(y[i-1]- y[i]) + signum(y[i+1]-y[i]);
        }
    }
    
    // now initialize the other n-1; actually not necessary because of default constructor
    for(int i=n; i<2*n-1; i++)
    {
        groupVec[i].active=false;
        groupVec[i].mergeLambda=-1;
        groupVec[i].mergeTo=-1;
    }
    
    // now calculate when the connected groups will hit and insert into 
    // multimap groupMove
    for(int i=0; i<n-1; ++i)
    {
        addConnection(i,i+1, 0);
    }
}


void FLSAClass::mergeGroups(int grpOne, int grpTwo, double lambda)
{
    // create a new group
    maxgroup++;
    groupDataNode one, two;
    
    one = groupVec[grpOne];
    two = groupVec[grpTwo];
    
    // activate the new group
    groupVec[maxgroup].active=true;
    groupVec[maxgroup].lambda=lambda;
    groupVec[maxgroup].mu = one.mu + (lambda-one.lambda) * one.deriv;
    groupVec[maxgroup].grpSize = one.grpSize + two.grpSize;
    groupVec[maxgroup].deriv = ((one.deriv * one.grpSize) + (two.deriv * two.grpSize))/groupVec[maxgroup].grpSize;
    
    // deactivate the two old groups
    deactivateGroup(grpOne, maxgroup, lambda);
    deactivateGroup(grpTwo, maxgroup, lambda);
    
    // find the neighbours of the new group and insert them
    vector<int> newNeighbour;
    newNeighbour = getNeighbours(grpOne, grpTwo);
    groupVec[maxgroup].neighbour.insert(groupVec[maxgroup].neighbour.begin(),newNeighbour.begin(), newNeighbour.end());
    newNeighbour = getNeighbours(grpTwo, grpOne);
    groupVec[maxgroup].neighbour.insert(groupVec[maxgroup].neighbour.begin(),newNeighbour.begin(), newNeighbour.end());
    
    // update the neighbours of the neighbours of the new group
    updateNeighbours(groupVec[maxgroup].neighbour, grpOne, maxgroup);
    updateNeighbours(groupVec[maxgroup].neighbour, grpTwo, maxgroup);

    // add the new connections
    groupDataNode newNode;
        
    newNode = groupVec[maxgroup];
    if(newNode.neighbour.size()>0)
    {
        for(unsigned int i=0; i< newNode.neighbour.size(); ++i)
        {
            addConnection(newNode.neighbour[i], maxgroup, lambda);
    
        }
    }
}



SEXP FLSAClass::prepSolTree(int numGrps)
{
    // prepare the list for the solution tree
    SEXP solTree;
    PROTECT(solTree=allocVector(VECSXP,5));
    
    // set the names of the components

    // SEXP names = getAttrib(solTree, R_NamesSymbol);  // Why bother if it is being clobbered next?
    SEXP names;
    PROTECT(names = allocVector(STRSXP,5));
    SET_STRING_ELT(names, 0, mkChar("mu"));
    SET_STRING_ELT(names, 1, mkChar("deriv"));
    SET_STRING_ELT(names, 2, mkChar("mergeLambda"));
    SET_STRING_ELT(names, 3, mkChar("mergeTo"));
    SET_STRING_ELT(names, 4, mkChar("numVars"));
    setAttrib(solTree, R_NamesSymbol, names);
    
    // now set the the class attribute to FLSA
    SEXP classStr;
    PROTECT(classStr = allocVector(STRSXP,1));
    SET_STRING_ELT(classStr, 0, mkChar("FLSA"));
    classgets(solTree, classStr);
    
    // generate teh vectors for the list of the right length
    SET_VECTOR_ELT(solTree,0, allocVector(REALSXP,numGrps));
    SET_VECTOR_ELT(solTree,1, allocVector(REALSXP,numGrps));
    SET_VECTOR_ELT(solTree,2, allocVector(REALSXP,numGrps));
    SET_VECTOR_ELT(solTree,3, allocVector(INTSXP,numGrps));
    SET_VECTOR_ELT(solTree,4, allocVector(INTSXP,1));
    
    UNPROTECT(3);
    return(solTree);

}


SEXP FLSAClass::solutionTree()
{
    SEXP solTree;
    // get the prepared tree
    PROTECT(solTree = prepSolTree(maxgroup+1)); 
    
    // fill in the data into the tree
    double* mu, *deriv, *mergeLambda;
    int* mergeTo, *numVars;
    mu = REAL(VECTOR_ELT(solTree,0));
    deriv = REAL(VECTOR_ELT(solTree,1));
    mergeLambda = REAL(VECTOR_ELT(solTree,2));
    mergeTo = INTEGER(VECTOR_ELT(solTree,3));
    for(int i=0; i<=maxgroup; ++i)
    {
        // set the 4 components
        mu[i] = groupVec[i].mu;
        deriv[i] = groupVec[i].deriv;
        mergeLambda[i] = groupVec[i].mergeLambda;
        mergeTo[i] = groupVec[i].mergeTo;
    }
    
    // save the number of variables
    numVars = INTEGER(VECTOR_ELT(solTree,4));
    numVars[0]=numVariables;
    
    // unprotect and return
    UNPROTECT(1);
    return(solTree);
}

/*
void FLSAClass::printGroupVec()
{
    for(int i=0; i<= maxgroup; ++i)
    {
        Rprintf("Mu[%d]: %f \n", i,groupVec[i].mu);
        Rprintf("Deriv[%d]: %f \n", i,groupVec[i].deriv);
        Rprintf("Size[%d]: %d\n", i, groupVec[i].grpSize);
    }
}

void FLSAClass::printGroupMove()
{
    typedef multimap<double, Connection>::const_iterator MI;
    int counter=0;
    Rprintf("Size: %d\n", groupMove.size());
    for(MI i=groupMove.begin(); i!= groupMove.end(); ++i)
    {
        counter++;
        Rprintf("%d: Lambda: %f; Connection (%d,%d) \n", counter,i->first, (i->second).first, (i->second).second);
    }
}
*/

/***********************************************
***
*** write a function that checks if everything with
*** y is ok
***
************************************************/

void FLSAClass::checkInput(SEXP y)
{
   // first check that y is a numeric vector
    if(!isNumeric(y))
    {
        error("y has to be a numeric vector");
    };
    
    // check that y is of at least length 2
    int len = LENGTH(y);
    if(len<2)
    {
        error("y has to be of length at least 2");
    };
}


extern "C"
{

SEXP FLSA(SEXP y)
{
    // initialize the object
    FLSAClass FLSAobj(y);
    pair<double, Connection> nextConn;
    
    // run until all groups have been fused (break when this is the case)
    while(true)
    {
//        FLSAobj.printGroupVec();
//        FLSAobj.printGroupMove();
        nextConn = FLSAobj.getNextConnection();
        // test if no more connections available
        if(nextConn.first==-1)
        {
            break;
        }
        else
        {
            FLSAobj.mergeGroups(nextConn.second.first, nextConn.second.second, nextConn.first);
        }
    }
    

    SEXP res;
    res = FLSAobj.solutionTree();
    
    return(res);

}


/******************************************************
***
*** function that given a solution tree and a vector with lambdas 
*** returns a matrix with the solution
*** lambdaVec is required to be sorted increasing
***
******************************************************/


SEXP FLSAexplicitSolution(SEXP solTree, SEXP lambdaR)
{
    SEXP resMat; // matrix for the result
    
    // data into the tree
    double* mu, *deriv, *mergeLambda;
    int* mergeTo, numVars;
    
    mu = REAL(VECTOR_ELT(solTree,0));
    deriv = REAL(VECTOR_ELT(solTree,1));
    mergeLambda = REAL(VECTOR_ELT(solTree,2));
    mergeTo = INTEGER(VECTOR_ELT(solTree,3));
    numVars = INTEGER(VECTOR_ELT(solTree,4))[0];
    
    // get a matrix in which to save the results 
    int lambdaLen = LENGTH(lambdaR);
    double* lambda = REAL(lambdaR);
    PROTECT(resMat = allocMatrix(REALSXP, lambdaLen, numVars));
    double* resMatVec = REAL(resMat); // for easier access
    int currGrp, currMatPos; // save the group that is currently in use
    double currLambda; // saves what the lambda level is 
    
    currMatPos = 0; //start at the beginning of the vector
    for(int i=0; i< numVars; ++i)
    {
        currGrp=i;
        currLambda=0;
        for(int j=0; j< lambdaLen; ++j)
        {
            while((mergeLambda[currGrp]< lambda[j]) && (mergeLambda[currGrp]!=-1)) // jump up to next group
            {
                currLambda = mergeLambda[currGrp];
                currGrp = mergeTo[currGrp];
            }
            // now am in the relevant group
            resMatVec[currMatPos] = mu[currGrp] + (lambda[j]-currLambda)*deriv[currGrp];
            currMatPos++;
        }
    }
    
    // unprotect and return
    UNPROTECT(1);
    return(resMat); // naming of rows and columns can be performed in R itself, easier there
}


}


