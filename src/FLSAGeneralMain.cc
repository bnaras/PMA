#include <iostream>
#include <fstream>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "PenaltyGraph.h"
#include "MaxFlowGraph.h"
#include "Groups.h"
#include <set>
#include "FLSAGeneral.h"

using namespace std;

// find the maximum of a vector of integers
int maxRIntVec(SEXP x)
{
//    Rprintf("Started max\n");
    int max=0;
    int len = LENGTH(x);
    int *xVec = INTEGER(x);
    for(int i =0; i< len; ++i)
    {
        if(xVec[i]>max)
        {
            max=xVec[i];
        }
    }
    return(max);
}

double maxRDoubleVec(SEXP x)
{
//    Rprintf("Started max\n");
    double max=-infinite;
    int len = LENGTH(x);
    double *xVec = REAL(x);
    for(int i =0; i< len; ++i)
    {
        if(xVec[i]>max)
        {
            max=xVec[i];
        }
    }
    return(max);
}


list<int> pointConn(int r, int c, int dimRow, int dimCol, int counter)
{
    list<int> foo;
    if(c>0) // connection to the left
    {
        foo.push_front(counter-dimRow);
    }
    if(c<(dimCol-1))// connection to the right
    {
        foo.push_front(counter+dimRow);
    }
    if(r>0) // connection above
    {
        foo.push_front(counter-1);
    }
    if(r<(dimRow-1))
    {
        foo.push_front(counter+1);
    }
    return(foo);
}

vector<list<int> > conn1Dim(int numOfNodes)
{
    vector<list<int> > conn(numOfNodes);
    list<int> foo;
    for(int i=0; i<numOfNodes; ++i)
    {
        foo.clear();
        if(i==0)
        {
            foo.push_front(1);
        }
        else if(i==numOfNodes-1)
        {
            foo.push_front(numOfNodes-2);
        }
        else
        {
            foo.push_front(i-1);
            foo.push_front(i+1);
        }
        conn[i]=foo;
    }
    return(conn);
}

vector<int> makeNodeNumbers1Dim(int numOfNodes)
{
    vector<int> nodeNumbers(numOfNodes);
    for(int i=0; i<numOfNodes; ++i)
    {
        nodeNumbers[i]=i;
    }
    return(nodeNumbers);
}

vector<double> readY(char* fileName)
{
    ifstream inputFile;
    list<double> help;
    double foo;
    int i=0;
    
    inputFile.open(fileName,ios::in);
    
    while(inputFile >> foo)
    {
        help.push_back(foo);
        ++i;
    }
    
    vector<double> y(help.size());
    for(unsigned int i=0; i<y.size(); ++i)
    {
        y[i]=help.front();
        help.pop_front();
    }
    inputFile.close();
    return(y);
}



extern "C" {

// main function that starts all the related FLSAGeneral classes

SEXP FLSAGeneralMain(SEXP connList, SEXP startValues, SEXP lambdas, SEXP maxSplitSize, SEXP verbose, SEXP thr, SEXP maxGrpNum)
{
    // find the highest nodeNumber
    int highNode = maxRIntVec(VECTOR_ELT(connList,0));
    SEXP sol;
    double highestLambda=infinite;
    
    if(IS_NUMERIC(lambdas))
    {
        highestLambda = maxRDoubleVec(lambdas);
    }
    
    FLSAGeneral FLSAGeneralObj(highNode, connList, startValues, maxSplitSize, verbose, thr, maxGrpNum, highestLambda);

    if(!IS_NUMERIC(lambdas))
    {
        sol=FLSAGeneralObj.solutionGraph();
    }
    else
    {
        sol =FLSAGeneralObj.solution(VECTOR_ELT(connList,0), lambdas);
    }
    return(sol);
}


SEXP FLSAGeneralExplicitSolution(SEXP solObj, SEXP nodes, SEXP lambdas)
{
    Groups groups(solObj);
    return(groups.solution(nodes, lambdas));
}


// generates the list of connections for the 2-dimensional fused lasso
// done in C++ to increase speed
// all variables are assumed checked for right format
// dimensions is a vector of two integers
// the nodes will be numbered starting at 0 along the columns (like an R matrix)
SEXP conn2Dim(SEXP dimensions)
{
    SEXP conn, bar;
    int dimRow = INTEGER(dimensions)[0];
    int dimCol = INTEGER(dimensions)[1];
    PROTECT(conn = allocVector(VECSXP, dimRow*dimCol));
    
    list<int> foo; // used to temporarily save the connections
    int counter = 0;
    // go through all gridpoints
    for(int j=0; j<dimCol; ++j)
    {
        for(int i=0; i<dimRow; ++i)
        {
            // get the connections of the current point
            foo=pointConn(i,j,dimRow,dimCol,counter);
            
            // now copy it into an R vector
            PROTECT(bar= allocVector(INTSXP, foo.size()));
            for(int k=0; k<LENGTH(bar); ++k)
            {
                INTEGER(bar)[k] = foo.front();
                foo.pop_front();
            }
            SET_VECTOR_ELT(conn,counter,bar);
            UNPROTECT(1);
            ++counter;
        }
    }
    
    UNPROTECT(1);
    return(conn);
}


} // end of extern C

