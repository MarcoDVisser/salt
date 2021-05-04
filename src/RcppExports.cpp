// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// resHyp
NumericVector resHyp(NumericVector wl, NumericVector R, NumericVector mu, NumericVector sigma);
RcppExport SEXP _salt_resHyp(SEXP wlSEXP, SEXP RSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type wl(wlSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R(RSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(resHyp(wl, R, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_salt_resHyp", (DL_FUNC) &_salt_resHyp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_salt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
