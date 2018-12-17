// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GxEScanC
int GxEScanC(Rcpp::IntegerVector& phenotype, Rcpp::NumericMatrix& covariates, Rcpp::IntegerVector& index, Rcpp::List& geneticData, std::string& outFilename, std::string& skippedFilename, int minmaf);
RcppExport SEXP _GxEScanR_GxEScanC(SEXP phenotypeSEXP, SEXP covariatesSEXP, SEXP indexSEXP, SEXP geneticDataSEXP, SEXP outFilenameSEXP, SEXP skippedFilenameSEXP, SEXP minmafSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type index(indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type geneticData(geneticDataSEXP);
    Rcpp::traits::input_parameter< std::string& >::type outFilename(outFilenameSEXP);
    Rcpp::traits::input_parameter< std::string& >::type skippedFilename(skippedFilenameSEXP);
    Rcpp::traits::input_parameter< int >::type minmaf(minmafSEXP);
    rcpp_result_gen = Rcpp::wrap(GxEScanC(phenotype, covariates, index, geneticData, outFilename, skippedFilename, minmaf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GxEScanR_GxEScanC", (DL_FUNC) &_GxEScanR_GxEScanC, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_GxEScanR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
