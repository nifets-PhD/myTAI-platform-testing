// Generated by cpp11: do not edit by hand
// clang-format off

#include <cpp11/R.hpp>
#include <Rcpp.h>
using namespace Rcpp;
#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// code.cpp
void fun();
extern "C" SEXP _myTAI_fun() {
  BEGIN_CPP11
    fun();
    return R_NilValue;
  END_CPP11
}

extern "C" {
/* .Call calls */
extern SEXP _myTAI_cpp_TAI(SEXP, SEXP);
extern SEXP _myTAI_cpp_bootMatrix(SEXP, SEXP, SEXP);
extern SEXP _myTAI_cpp_geom_mean(SEXP);
extern SEXP _myTAI_cpp_harmonic_mean(SEXP);
extern SEXP _myTAI_cpp_omitMatrix(SEXP, SEXP);
extern SEXP _myTAI_cpp_pMatrix(SEXP, SEXP);
extern SEXP _myTAI_cpp_std_error(SEXP);
extern SEXP _myTAI_rcpp_boottei_parallel(SEXP, SEXP, SEXP, SEXP);
extern SEXP _myTAI_rcpp_pMatrix_parallel(SEXP, SEXP, SEXP);
extern SEXP _myTAI_rcpp_pStrata_parallel(SEXP, SEXP, SEXP, SEXP);
extern SEXP _myTAI_rcpp_tei_parallel(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_myTAI_cpp_TAI",               (DL_FUNC) &_myTAI_cpp_TAI,               2},
    {"_myTAI_cpp_bootMatrix",        (DL_FUNC) &_myTAI_cpp_bootMatrix,        3},
    {"_myTAI_cpp_geom_mean",         (DL_FUNC) &_myTAI_cpp_geom_mean,         1},
    {"_myTAI_cpp_harmonic_mean",     (DL_FUNC) &_myTAI_cpp_harmonic_mean,     1},
    {"_myTAI_cpp_omitMatrix",        (DL_FUNC) &_myTAI_cpp_omitMatrix,        2},
    {"_myTAI_cpp_pMatrix",           (DL_FUNC) &_myTAI_cpp_pMatrix,           2},
    {"_myTAI_cpp_std_error",         (DL_FUNC) &_myTAI_cpp_std_error,         1},
    {"_myTAI_fun",                   (DL_FUNC) &_myTAI_fun,                   0},
    {"_myTAI_rcpp_boottei_parallel", (DL_FUNC) &_myTAI_rcpp_boottei_parallel, 4},
    {"_myTAI_rcpp_pMatrix_parallel", (DL_FUNC) &_myTAI_rcpp_pMatrix_parallel, 3},
    {"_myTAI_rcpp_pStrata_parallel", (DL_FUNC) &_myTAI_rcpp_pStrata_parallel, 4},
    {"_myTAI_rcpp_tei_parallel",     (DL_FUNC) &_myTAI_rcpp_tei_parallel,     3},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_myTAI(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
