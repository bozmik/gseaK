// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ES
NumericMatrix ES(NumericVector stat, CharacterVector geneOrd, double Pmiss, double Nr, CharacterVector gs);
RcppExport SEXP gseak_ES(SEXP statSEXP, SEXP geneOrdSEXP, SEXP PmissSEXP, SEXP NrSEXP, SEXP gsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type geneOrd(geneOrdSEXP);
    Rcpp::traits::input_parameter< double >::type Pmiss(PmissSEXP);
    Rcpp::traits::input_parameter< double >::type Nr(NrSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type gs(gsSEXP);
    rcpp_result_gen = Rcpp::wrap(ES(stat, geneOrd, Pmiss, Nr, gs));
    return rcpp_result_gen;
END_RCPP
}
// ginGS
List ginGS(CharacterVector g, CharacterVector gs, NumericVector stat);
RcppExport SEXP gseak_ginGS(SEXP gSEXP, SEXP gsSEXP, SEXP statSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type gs(gsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    rcpp_result_gen = Rcpp::wrap(ginGS(g, gs, stat));
    return rcpp_result_gen;
END_RCPP
}
// permutation
NumericMatrix permutation(NumericVector x, int n);
RcppExport SEXP gseak_permutation(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(permutation(x, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"gseak_ES", (DL_FUNC) &gseak_ES, 5},
    {"gseak_ginGS", (DL_FUNC) &gseak_ginGS, 3},
    {"gseak_permutation", (DL_FUNC) &gseak_permutation, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_gseak(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
