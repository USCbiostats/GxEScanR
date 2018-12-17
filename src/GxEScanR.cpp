// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <iostream>
#include <fstream>
#include <RcppArmadillo.h>

// [[Rcpp::export]]
int GxEScanC(Rcpp::IntegerVector &phenotype, Rcpp::NumericMatrix &covariates, Rcpp::IntegerVector &index,
             Rcpp::List &geneticData, std::string &outFilename, std::string &skippedFilename, int minmaf) {
  return 0;
}