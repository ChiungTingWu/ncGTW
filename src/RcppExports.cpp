// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// graphCut
NumericVector graphCut(NumericMatrix termW, NumericMatrix edges);
RcppExport SEXP _ncGTW_graphCut(SEXP termWSEXP, SEXP edgesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type termW(termWSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type edges(edgesSEXP);
    rcpp_result_gen = Rcpp::wrap(graphCut(termW, edges));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ncGTW_graphCut", (DL_FUNC) &_ncGTW_graphCut, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ncGTW(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
