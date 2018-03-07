// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <iostream>
#include <fstream>
#include <Rcpp.h>
#include "Subject.h"
#include "GeneticData.h"
#include "BinaryDosage.h"
#include "MatrixFunctions.h"
#include "GxEDataset.h"

void ConvertCovariates(Rcpp::List subjectData, std::vector<double> &cov, int numSubjects, int numCov) {
  std::vector<double> covR  = Rcpp::as<std::vector<double> >(subjectData["covariates"]);

  cov.resize(covR.size());
  Rcpp::Rcout << numSubjects << '\t' << numCov << '\t' << covR.size() << std::endl;
  Transpose(&covR[0], &cov[0], (unsigned int)numSubjects, (unsigned int)numCov);
}
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
  std::vector<std::vector<double> > probabilities;
  std::vector<int> subjectOrder = Rcpp::as<std::vector<int> >(subjectData["gLoc"]);
  std::vector<double> phenotype = Rcpp::as<std::vector<double> >(subjectData["phenotypes"]);
  std::vector<bool> missingPhenotype;
  Rcpp::List subList = geneticInfo["subjects"];
  Rcpp::DataFrame subInfo = Rcpp::as<Rcpp::DataFrame>(subList["Info"]);
  std::vector<std::string> subjectID = Rcpp::as<std::vector<std::string> >(subInfo["IID"]);
  std::vector<double> covariates;
  std::vector<bool> missingCov;
  std::vector<bool> filter;
  std::vector<double> geneticValues;
  int numSubjectsUsed, numCov;
  int i;
  double *d, *p0, *p1, *p2;
  Rcpp::List res;
  
  numSubjectsUsed = subjectOrder.size();
  numCov = Rcpp::as<std::vector<double> >(subjectData["covariates"]).size() / numSubjectsUsed;
  ConvertCovariates(subjectData, covariates, numCov, numSubjectsUsed);
  Rcpp::Rcout << "Number of covariate values\t" << numCov << std::endl;
  format = (int)geneticInfo["format"];
  subversion = (int)geneticInfo["version"];
  numSubjects = (int)geneticInfo["numSubjects"];
  numSNPs = (int)geneticInfo["numSNPs"];
  gFilename = Rcpp::as<std::string>(geneticInfo["geneticFile"]);

  if (subversion == 1)
    geneticValues.resize(numSubjectsUsed);
  else
    geneticValues.resize(4 * numSubjectsUsed);
  
  missingCov.assign(covariates.size(), false);
  missingPhenotype.assign(phenotype.size(), false);
  filter.assign(numSubjectsUsed, true);
  
  if (format == 1 && subversion == 1)
    geneticData = new CBinaryDosageFormat1_1(gFilename, numSubjects, numSNPs);
  else if (format == 1 && subversion == 2)
    geneticData = new CBinaryDosageFormat1_2(gFilename, numSubjects, numSNPs);
  else if (format == 2 && subversion == 1)
    geneticData = new CBinaryDosageFormat2_1(gFilename, numSubjects, numSNPs);
  else if (format == 2 && subversion == 2)
    geneticData = new CBinaryDosageFormat2_2(gFilename, numSubjects, numSNPs);
  else if (format == 3 && subversion == 1)
    geneticData = new CBinaryDosageFormat3_1(gFilename, numSubjects, numSNPs);
  else if (format == 3 && subversion == 2)
    geneticData = new CBinaryDosageFormat3_2(gFilename, numSubjects, numSNPs);
  else if (format == 4 && subversion == 2)
    geneticData = new CBinaryDosageFormat4_2(gFilename, numSubjects, numSNPs);
  
  if (geneticData != NULL) {
    if (geneticData->GetFirst()) {
      Rcpp::Rcout << "GetFirst failure" << std::endl;
      Rcpp::Rcout << geneticData->ErrorMessage() << std::endl;
    } else {
      dosages = geneticData->Dosages();
      d = &geneticValues[0];
      if (subversion == 2) {
        p0 = d + numSubjectsUsed;
        p1 = p0 + numSubjectsUsed;
        p2 = p1 + numSubjectsUsed;
        for (i = 0; i < numSubjectsUsed; ++i, ++d, ++p0, ++p1, ++p2) {
          *d = geneticData->Dosages()[subjectOrder[i] - 1];
          *p0 = geneticData->Probabilities()[0][subjectOrder[i] - 1];
          *p1 = geneticData->Probabilities()[1][subjectOrder[i] - 1];
          *p2 = geneticData->Probabilities()[2][subjectOrder[i] - 1];
        }
      } else {
        for (i = 0; i < numSubjectsUsed; ++i, ++d)
          *d = geneticData->Dosages()[subjectOrder[i] - 1];
      }
      probabilities = geneticData->Probabilities();
    }

    if (geneticData->GetNext()) {
      Rcpp::Rcout << "GetNext failure" << std::endl;
      Rcpp::Rcout << geneticData->ErrorMessage() << std::endl;
    } else {
      dosages = geneticData->Dosages();
      d = &geneticValues[0];
      if (subversion == 2) {
        p0 = d + numSubjectsUsed;
        p1 = p0 + numSubjectsUsed;
        p2 = p1 + numSubjectsUsed;
        for (i = 0; i < numSubjectsUsed; ++i, ++d, ++p0, ++p1, ++p2) {
          *d = geneticData->Dosages()[subjectOrder[i] - 1];
          *p0 = geneticData->Probabilities()[0][subjectOrder[i] - 1];
          *p1 = geneticData->Probabilities()[1][subjectOrder[i] - 1];
          *p2 = geneticData->Probabilities()[2][subjectOrder[i] - 1];
        }
      } else {
        for (i = 0; i < numSubjectsUsed; ++i, ++d)
          *d = geneticData->Dosages()[subjectOrder[i] - 1];
      }
      probabilities = geneticData->Probabilities();
    }
  }
  if (geneticData)
    delete geneticData;
  
  res = Rcpp::List::create(Rcpp::Named("Dosages") = dosages,
                           Rcpp::Named("GeneValues") = geneticValues,
                           Rcpp::Named("Phenotypes") = phenotype,
                           Rcpp::Named("Covariates") = covariates,
                           Rcpp::Named("Probabilities") = probabilities);
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