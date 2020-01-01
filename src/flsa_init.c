#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call entry points */
extern SEXP FLSA(SEXP y); 
extern SEXP FLSAexplicitSolution(SEXP solTree, SEXP lambdaVec); 
extern SEXP FLSAGeneralMain(SEXP connList, SEXP startValues, SEXP lambdas, SEXP maxSplitSize, SEXP verbose, SEXP thr, SEXP maxGrpNum);
extern SEXP FLSAGeneralExplicitSolution(SEXP solObj, SEXP nodes, SEXP lambdas);
extern SEXP conn2Dim(SEXP dimensions);

static const R_CallMethodDef CallEntries[] = {
					      {"FLSA", (DL_FUNC) &FLSA, 1},
					      {"FLSAexplicitSolution", (DL_FUNC) &FLSAexplicitSolution, 2},
					      {"FLSAGeneralMain", (DL_FUNC) &FLSAGeneralMain, 7},
					      {"FLSAGeneralExplicitSolution", (DL_FUNC) &FLSAGeneralExplicitSolution, 3},
					      {"conn2Dim", (DL_FUNC) &conn2Dim, 1},
					      {NULL, NULL, 0}
};
    
void R_init_flsa(DllInfo *info) {
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
  

