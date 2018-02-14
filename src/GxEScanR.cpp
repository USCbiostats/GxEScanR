// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <iostream>
#include <fstream>
#include <Rcpp.h>
#include "Subject.h"
#include "GeneticData.h"
#include "BinaryDosage.h"

//' Function to fit models scanning over genotypes
//' 
//' Function to fit selected models over genotypes. Results from
//' these models are then processed by GxETest to display the results
//' for these one step tests along with two step tests derived from
//' the one step tests.
//' 
//' @param subjectData
//' List returned from SubsetSubjects
//' @param geneticInfo
//' List returned from one of the functions to get the required information
//' about the source of genetic data.
//' @return
//' 0 success
//' 1 failure
//' @export
// [[Rcpp::export]]
Rcpp::List GxEScanC(Rcpp::List subjectData, Rcpp::List geneticInfo) {
  CGeneticData *geneticData = NULL;
  int format, subversion;
  int numSubjects, numSNPs;
  std::string gFilename;
  std::vector<double> dosages;
  Rcpp::List res;
  
  format = (int)geneticInfo["format"];
  subversion = (int)geneticInfo["version"];
  Rcpp::Rcout << "Format:\t"<< format << '.' << subversion << std::endl;
  numSubjects = (int)geneticInfo["numSubjects"];
  numSNPs = (int)geneticInfo["numSNPs"];
  gFilename = Rcpp::as<std::string>(geneticInfo["geneticFile"]);
  Rcpp::Rcout << "Filename:\t" << gFilename << std::endl;
  Rcpp::Rcout << "Number of subjects:\t" << numSubjects << std::endl;
  Rcpp::Rcout << "NUmber of SNPs:\t" << numSNPs << std::endl; 
  if (format == 1 && subversion == 1) {
    geneticData = new CBinaryDosageFormat1_1(gFilename, numSubjects, numSNPs);
    if (geneticData->GetFirst()) {
      Rcpp::Rcout << "GetFirst failure" << std::endl;
      Rcpp::Rcout << geneticData->ErrorMessage() << std::endl;
    }
    dosages = geneticData->Dosages();
    
  }
  if (geneticData)
    delete geneticData;
  
  res = Rcpp::List::create(Rcpp::Named("Dosages") = dosages);
  return res;
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