#include <R.h>
#include <Rinternals.h>
#include <cstdlib> // for NULL
#include <R_ext/Rdynload.h>


extern "C"
{
    // declaring functions for exporting
    SEXP FLSA(SEXP y); 
    SEXP FLSAexplicitSolution(SEXP solTree, SEXP lambdaVec); 
    SEXP FLSAGeneralMain(SEXP connList, SEXP startValues, SEXP lambdas, SEXP maxSplitSize, SEXP thr, SEXP maxGrpNum);
    SEXP FLSAGeneralExplicitSolution(SEXP solObj, SEXP nodes, SEXP lambdas);
    SEXP conn2Dim(SEXP dimensions);
    
    void R_init_flsa(DllInfo *info)
    {
        R_CallMethodDef  callMethods[] = {
            {"FLSA", (void* (*)()) &FLSA, 1},
            {"FLSAexplicitSolution", (void* (*)()) &FLSAexplicitSolution,2},
            {"FLSAGeneralMain", (void* (*)()) &FLSAGeneralMain, 7},
            {"FLSAGeneralExplicitSolution", (void* (*)()) &FLSAGeneralExplicitSolution, 3},
            {"conn2Dim", (void* (*)()) &conn2Dim, 1},
            {NULL, NULL, 0}
        };
        
        R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
    }

    void R_unload_flsa(DllInfo *info)
    {
    
    }
}


