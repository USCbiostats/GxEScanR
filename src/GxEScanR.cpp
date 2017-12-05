// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "Subject.h"
#include "GeneticData.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

//' Function to fit models scanning over genotypes
//' 
//' Function to fit selected models over genotypes. Results from
//' these models are then processed by GxETest to display the results
//' for these one step tests along with two step tests derived from
//' the one step tests.
//' 
//' @param subjectData
//' Data frame with family and subject IDs in first two columns
//' phenotype in third columna
//' and covariates in remaining columns
//' @return
//' 0 success
//' 1 failure
//' @export
// [[Rcpp::export]]
int GxEScan(Rcpp::DataFrame subjectData, Rcpp::List geneticData) {
//  Rcpp::Rcout << "Checking subject data" << std::endl;
  Rcpp::Rcout << "Checking subjectData:\t\t";
  if (TestSubjectData(subjectData)) {
    Rcpp::Rcerr << "Failed" << std::endl;
    return 1;
  }
  Rcpp::Rcout << "Passed" << std::endl;
  Rcpp::Rcout << "Checking genetic data:\t\t";
  if (TestGeneticData(geneticData)) {
    Rcpp::Rcerr << "Failed" << std::endl;
    return 1;
  }
  Rcpp::Rcout << "Passed" << std::endl;
  Rcpp::Rcout << "Checking model selection" << std::endl;
  Rcpp::Rcout << "Fitting models" << std::endl;
  Rcpp::Rcout << "Complete" << std::endl;
  return 0;
}

//' Function to display the results from the scans performed by GxEScan
//' 
//' Function to display the results from the scans performed by GxEScan.
//' The results include Manhattan and qq-plot along with tables of the
//' top hits for each test. The results are saved to an Excel spreadsheet.
//' 
//' @return
//' 0 success
//' 1 failure
//' @export
// [[Rcpp::export]]
int GxETest() {
  Rcpp::Rcout << "Checking test selection" << std::endl;
  Rcpp::Rcout << "Writing summary page" << std::endl;
  Rcpp::Rcout << "Processing one step tests" << std::endl;
  Rcpp::Rcout << "Processing two step tests" << std::endl;
  Rcpp::Rcout << "Complete" << std::endl;
  return 0;
}