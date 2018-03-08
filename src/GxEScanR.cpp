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
  Transpose(covR.data(), cov.data(), (unsigned int)numSubjects, (unsigned int)numCov);
}

void WriteResults(std::ostream &outfile, const int maxRet, const CGxEPolytomousDataset &gxeData) {
  int betaLoc, nCov;
  
  if (maxRet & 0x7FF)
    return;
  
  nCov = gxeData.NumCovariates();
  if (maxRet & 0x01) {
    outfile << "NA\tNA";
  } else {
    betaLoc = nCov - 1;
    outfile << gxeData.BetaD_GE()[betaLoc] / gxeData.Mean()[betaLoc] << '\t'
            << gxeData.BetaD_GE()[betaLoc] * sqrt(gxeData.InverseInformationD_GE()[nCov * nCov - 1]);
  }
  outfile << std::endl;
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
//' @param outputFilename
//' Name of output file
//' @return
//' 0 success
//' 1 failure
//' @export
// [[Rcpp::export]]
Rcpp::List GxEScanC(Rcpp::List subjectData, Rcpp::List geneticInfo, std::string outputFilename) {
  CGeneticData *geneticData = NULL;
  int format, subversion;
  int numSubjects, numSNPs;
  std::string gFilename;
  std::vector<double> dosages;
  std::vector<std::vector<double> > probabilities;
  std::vector<int> subjectOrder = Rcpp::as<std::vector<int> >(subjectData["gLoc"]);
  std::vector<double> phenotype = Rcpp::as<std::vector<double> >(subjectData["phenotypes"]);
  std::vector<char> caseControl;
  std::vector<char> missingPhenotype;
  Rcpp::List subList = geneticInfo["subjects"];
  Rcpp::DataFrame subInfo = Rcpp::as<Rcpp::DataFrame>(subList["Info"]);
  std::vector<std::string> subjectID = Rcpp::as<std::vector<std::string> >(subInfo["IID"]);
  std::vector<double> covariates;
  std::vector<char> missingCov;
  std::vector<char> missingGene;
  std::vector<char> filter;
  std::vector<double> geneticValues;
  int numSubjectsUsed, numCov;
  int i, j;
  unsigned int ui;
  int retVal;
  double *d, *p0, *p1, *p2;
  std::ofstream outfile(outputFilename.c_str());
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
  
  missingPhenotype.assign(phenotype.size(), false);
  caseControl.assign(phenotype.size(), false);
  for (ui = 0; ui < phenotype.size(); ++ui) {
    if (phenotype[ui] == 1)
      caseControl[ui] = true;
    else if (phenotype[ui] != 0)
      missingPhenotype[ui] = true;
  }
  missingCov.assign(covariates.size(), false);
  missingGene.assign(numSubjectsUsed, false);
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
/*
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
*/
  }
  
  CGxEPolytomousDataset gxeData(numSubjectsUsed, caseControl.data(), (bool *)missingPhenotype.data(),
                                numCov, covariates.data(), (bool *)missingCov.data(), numCov, (bool *)filter.data());
  if (gxeData.Initialize() == false)
    Rcpp::Rcout << "Initialization failed - " << gxeData.ErrorString() << std::endl;
  if (subversion == 2)
    gxeData.AssignGene(geneticValues.data(), (bool *)missingGene.data(), true, true);
  else
    gxeData.AssignGene(geneticValues.data(), (bool *)missingGene.data(), true, false);
  gxeData.UpdateGene();
  retVal = gxeData.FitModels();
  WriteResults(outfile, retVal, gxeData);
  
  if (geneticData)
    delete geneticData;

  std::vector<double> covInt;
  covInt.assign(numSubjectsUsed * gxeData.NumParameters(), -9);
  if (gxeData.CompleteCovariates() != NULL) {
    for (i = 0; i < numSubjectsUsed * gxeData.NumParameters(); ++i)
      covInt[i] = gxeData.CompleteCovariates()[i];
  }
  /*
  std::ofstream outfile("GxEDataset.txt");
  for (i = 0; i < numSubjects; ++i) {
    outfile << phenotype[i];
    for (j = 0; j < gxeData.NumParameters(); ++j)
      outfile << '\t' << covInt[i * gxeData.NumParameters() + j];
    outfile << std::endl;
  }
  outfile.close();
 */
  Rcpp::Rcout << "Number of parameters\t" << gxeData.NumParameters() << std::endl;
  std::vector<double> covMean;
  covMean.resize(gxeData.NumParameters());
  for (i = 0; i < gxeData.NumParameters(); ++i)
    covMean[i] = gxeData.Mean()[i];
  std::vector<double> beta;
  beta.resize(gxeData.NumCovariates() + 2);
  for (ui = 0; ui < gxeData.NumCovariates() + 2; ++ui)
    beta[ui] = gxeData.BetaD_GxE()[ui];
  std::vector<double> invInfo;
  invInfo.resize((gxeData.NumCovariates() + 2) * (gxeData.NumCovariates() + 2));
  for (ui = 0; ui < (gxeData.NumCovariates() + 2) * (gxeData.NumCovariates() + 2); ++ui)
    invInfo[ui] = gxeData.InverseInformationD_GxE()[ui];
  res = Rcpp::List::create(Rcpp::Named("Dosages") = dosages,
                           Rcpp::Named("GeneValues") = geneticValues,
                           Rcpp::Named("Phenotypes") = phenotype,
                           Rcpp::Named("Covariates") = covariates,
                           Rcpp::Named("CovInt") = covInt,
                           Rcpp::Named("CovMean") = covMean,
                           Rcpp::Named("MaxRet") = retVal,
                           Rcpp::Named("Beta") = beta,
                           Rcpp::Named("InvInfo") = invInfo,
                           Rcpp::Named("Probabilities") = probabilities);
  outfile.close();
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