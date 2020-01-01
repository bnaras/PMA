#include <R.h>
#include <Rinternals.h>
#include <vector>
#include <map>

using namespace std;

typedef pair<int,int> Connection;


struct groupDataNode {
    bool active;
    double mu, lambda, deriv, mergeLambda;
    int grpSize, mergeTo;
    vector<int> neighbour;
};

// make a class for all the groups

class FLSAClass {
    // used to describe nodes that are connected

    vector<groupDataNode> groupVec;
    multimap<double, Connection> groupMove; // multimap to store the times when groups hit each other
    int maxgroup;
    int numVariables;
    
    void checkInput(SEXP y); // checks if y has the right input format
    void addConnection(int grpOne, int grpTwo, double lambda);
    double getCurrMu(groupDataNode x, double lambda) {return(x.mu + (lambda-x.lambda)*x.deriv);};
    vector<int> getNeighbours(int grpNum, int exclGrp); // get the neighbours of grpNum, excluding exclGrp
    void updateNeighbours(vector<int> updateGrp, int oldGrp, int newGrp); // changes the grpNumber of the neighbours from old to new
    void deactivateGroup(int grpNum, int newGrp, double lambda); // deactivate group grpNum, which merges into newGrp at lambda
    SEXP prepSolTree(int numGrps); // prepare the list for the solution tree
    
public:
    FLSAClass(SEXP y); // Constructor
    void mergeGroups(int grpOne, int grpTwo, double lambda); // merges two groups
    pair<double, Connection> getNextConnection(); // returns the connection with the smallest lambda value
    SEXP solutionTree(); // returns the found solution so far in the form of a tree described in vector format
//    void printGroupVec();
//    void printGroupMove();
};


extern "C"
{

//main function for FLSA
SEXP FLSA(SEXP y); 

// function that for a vector of lambdas explicitly calculates and returns the solution
SEXP FLSAexplicitSolution(SEXP solTree, SEXP lambdaVec); 
}
