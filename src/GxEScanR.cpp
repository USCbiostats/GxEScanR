// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <iostream>
#include <fstream>
#include <Rcpp.h>
#include "Subject.h"
#include "GeneticData.h"
#include "BinaryDosage.h"
#include "Impute2.h"
#include "MatrixFunctions.h"
#include "GxEDataset.h"

CGeneticData *OpenBinaryDosageFile(const Rcpp::List &geneticInfo) {
  CGeneticData *geneticData = NULL;
  int format, subversion, numSubjects, numSNPs;
  std::string filename;
  
  format = (int)geneticInfo["format"];
  subversion = (int)geneticInfo["version"];
  numSubjects = (int)geneticInfo["numSubjects"];
  numSNPs = (int)geneticInfo["numSNPs"];
  filename = Rcpp::as<std::string>(geneticInfo["filename"]);
  
  if (format == 1 && subversion == 1)
    geneticData = new CBinaryDosageFormat1_1(filename, numSubjects, numSNPs);
  else if (format == 1 && subversion == 2)
    geneticData = new CBinaryDosageFormat1_2(filename, numSubjects, numSNPs);
  else if (format == 2 && subversion == 1)
    geneticData = new CBinaryDosageFormat2_1(filename, numSubjects, numSNPs);
  else if (format == 2 && subversion == 2)
    geneticData = new CBinaryDosageFormat2_2(filename, numSubjects, numSNPs);
  else if (format == 3 && subversion == 1)
    geneticData = new CBinaryDosageFormat3_1(filename, numSubjects, numSNPs);
  else if (format == 3 && subversion == 2)
    geneticData = new CBinaryDosageFormat3_2(filename, numSubjects, numSNPs);
  else if (format == 4 && subversion == 1)
    geneticData = new CBinaryDosageFormat4_1(filename, numSubjects, numSNPs);
  else if (format == 4 && subversion == 2)
    geneticData = new CBinaryDosageFormat4_2(filename, numSubjects, numSNPs);
  else
    Rcpp::stop("Unknown binary dosage file format");
  
  return geneticData;
}

CGeneticData *OpenGeneticData(const Rcpp::List &geneticInfo) {
  std::string geneticFileType;

  geneticFileType = Rcpp::as<std::string>(geneticInfo["filetype"]);
  if (geneticFileType == "BinaryDosage")
    return OpenBinaryDosageFile(geneticInfo);
  if (geneticFileType == "Impute2")
    return OpenImpute2File(geneticInfo);
  Rcpp::stop("Unknown genetic data type");
  
  return NULL;
}

void ConvertCovariates(Rcpp::List subjectData, std::vector<double> &cov, int numSubjects, int numCov) {
  std::vector<double> covR  = Rcpp::as<std::vector<double> >(subjectData["covariates"]);

  cov.resize(covR.size());
  Transpose(covR.data(), cov.data(), (unsigned int)numSubjects, (unsigned int)numCov);
}

void WriteResults(std::ostream &outfile, const int maxRet, const CGxEPolytomousDataset &gxeData) {
  int betaLoc, nCov;
  double a, b, c, x, y;
  
  if ((maxRet & 0x7ff) == 0x7ff)
    return;
  
  nCov = gxeData.NumCovariates();
//  Rcpp::Rcout << "Number of Covariates:\t" << nCov << std::endl;
  if (maxRet & 0x01) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov;
    outfile << gxeData.BetaD_GE()[betaLoc] << '\t'
            << gxeData.BetaD_GE()[betaLoc] / sqrt(gxeData.InverseInformationD_GE()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x02) {
    outfile << "NA\tNA\tNA\t";
  } else {
    betaLoc = nCov + 1;
    a = gxeData.InverseInformationD_GxE()[(betaLoc + 2) * (betaLoc - 1)];
    b = gxeData.InverseInformationD_GxE()[(betaLoc + 2) * (betaLoc - 1) + 1];
    c = gxeData.InverseInformationD_GxE()[(betaLoc + 2) * betaLoc];
    x = gxeData.BetaD_GxE()[betaLoc - 1];
    y = gxeData.BetaD_GxE()[betaLoc];
    outfile << y / gxeData.Mean()[betaLoc - 2] << '\t'
            << y / sqrt(c) << '\t'
            << (c * x * x - (b + b) * x * y + a * y * y) / (a * c - b * b) << '\t';
  }
  if (maxRet &0x04) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov - 1;
    outfile << gxeData.BetaG_E()[betaLoc] / gxeData.Mean()[betaLoc] << '\t'
            << gxeData.BetaG_E()[betaLoc] / sqrt(gxeData.InverseInformationG_E()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x08) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov - 1;
    outfile << gxeData.BetaCaseOnly()[betaLoc] / gxeData.Mean()[betaLoc] << '\t'
            << gxeData.BetaCaseOnly()[betaLoc] / sqrt(gxeData.InverseInformationCO()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x10) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov - 1;
    outfile << gxeData.BetaCntlOnly()[betaLoc] / gxeData.Mean()[betaLoc] << '\t'
            << gxeData.BetaCntlOnly()[betaLoc] / sqrt(gxeData.InverseInformationCntlOnly()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x20) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov + nCov - 2;
    outfile << gxeData.BetaPolytomousG_E()[betaLoc] / gxeData.Mean()[nCov - 1] << '\t'
            << gxeData.BetaPolytomousG_E()[betaLoc] / sqrt(gxeData.InverseInformationPolytomousG_E()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x40) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov + nCov - 2;
    outfile << gxeData.BetaPolytomousCaseOnly()[betaLoc] / gxeData.Mean()[nCov - 1] << '\t'
            << gxeData.BetaPolytomousCaseOnly()[betaLoc] / sqrt(gxeData.InverseInformationPolytomousCaseOnly()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x80) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov + nCov - 2;
    outfile << gxeData.BetaPolytomousControlOnly()[betaLoc] / gxeData.Mean()[nCov - 1] << '\t'
            << gxeData.BetaPolytomousControlOnly()[betaLoc] / sqrt(gxeData.InverseInformationPolytomousControlOnly()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x100) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov;
    outfile << gxeData.BetaRestrictedPolytomousG_E()[betaLoc] / gxeData.Mean()[nCov - 1] << '\t'
            << gxeData.BetaRestrictedPolytomousG_E()[betaLoc] / sqrt(gxeData.InverseInformationRestrictedPolytomousG_E()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x200) {
    outfile << "NA\tNA\t";
  } else {
    betaLoc = nCov;
    outfile << gxeData.BetaRestrictedPolytomousCaseOnly()[betaLoc] / gxeData.Mean()[nCov - 1] << '\t'
            << gxeData.BetaRestrictedPolytomousCaseOnly()[betaLoc] / sqrt(gxeData.InverseInformationRestrictedPolytomousCaseOnly()[(betaLoc + 2) * betaLoc]) << '\t';
  }
  if (maxRet &0x400) {
    outfile << "NA\tNA";
  } else {
    betaLoc = nCov;
    outfile << gxeData.BetaRestrictedPolytomousControlOnly()[betaLoc] / gxeData.Mean()[nCov - 1] << '\t'
            << gxeData.BetaRestrictedPolytomousControlOnly()[betaLoc] / sqrt(gxeData.InverseInformationRestrictedPolytomousControlOnly()[(betaLoc + 2) * betaLoc]);
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
int GxEScanC(Rcpp::List subjectData, Rcpp::List geneticInfo, std::string outputFilename) {
  CGeneticData *geneticData = NULL;
  std::vector<double> phenotype = Rcpp::as<std::vector<double> >(subjectData["phenotype"]);
  std::vector<double> covariates;
  std::vector<int> geneIndex = Rcpp::as<std::vector<int> >(subjectData["geneIndex"]);
  int numSubjectsUsed, numCov;
  // Data needed for GxEDataset
  std::vector<char> missingPhenotype;
  std::vector<char> caseControl; // Already have variable named phenotype
  std::vector<char> missingCov;
  std::vector<char> missingGene;
  std::vector<char> filter;
  std::vector<double> geneticValues;
  double *d, *p0, *p1, *p2;
  
  std::ofstream outfile(outputFilename.c_str());
  int retVal;

  int i, j; 
  unsigned int ui;
  /*
  int format, subversion;
  int numSubjects, numSNPs;
  std::string geneticFileType;
  std::string gFilename;
  std::vector<double> dosages;
  std::vector<std::vector<double> > probabilities;
  Rcpp::List subList = geneticInfo["subjects"];
  Rcpp::DataFrame subInfo = Rcpp::as<Rcpp::DataFrame>(subList["Info"]);
  std::vector<std::string> subjectID = Rcpp::as<std::vector<std::string> >(subInfo["IID"]);
  std::vector<double> covariates;
  std::vector<double> geneticValues;
  int numSubjectsUsed, numCov;
  int i, j;
  unsigned int ui;
  int retVal;
  double *d, *p0, *p1, *p2;
  std::ofstream outfile(outputFilename.c_str());
  Rcpp::List res;
*/ 
  if (!outfile.good()) {
    Rcpp::stop("Unable to open output file");
  }
  
  Rcpp::DataFrame snpData = Rcpp::as<Rcpp::DataFrame>(geneticInfo["snps"]);
  std::vector<std::string> snpName = Rcpp::as<std::vector<std::string> >(snpData["SNP"]);
  std::vector<std::string> chrName = Rcpp::as<std::vector<std::string> >(snpData["CHR"]);
  std::vector<std::string> a1Name = Rcpp::as<std::vector<std::string> >(snpData["A1"]);
  std::vector<std::string> a2Name = Rcpp::as<std::vector<std::string> >(snpData["A2"]);
  std::vector<int> bp = Rcpp::as<std::vector<int> >(snpData["BP"]);
  if (chrName[0] != "")
    outfile << "CHR\t";
  if (snpName[0] != "")
    outfile << "SNP\t";
  if (bp[0] != 0)
    outfile << "BP\t";
  if (a1Name[0] != "")
    outfile << "A1\tA2\t";
//  Rcpp::Rcout << snpName[0] << '\t' << chrName[0] << '\t' << a1Name[0] << '\t' << a2Name[0] << '\t' << bp[0] << std::endl;
  outfile << "Cases\tControls\tBetaG\tzG\tBetaGxE\tzGxE\tChi2df\t";
  outfile << "Beta_HWGE\tz_HWGE\tBeta_HWCase\tz_HWCase\tBeta_HWCtrl\tz_HWCtrl\t";
  outfile << "Beta_PolyGE\tz_PolyGE\tBeta_PolyCase\tz_PolyCase\tBeta_PolyCtrl\tzPolyCtrl\t";
  outfile << "Beta_RPolyGE\tz_RPolyGE\tBeta_RPolyCase\tz_RPolyCase\tBeta_RPolyCtrl\tz_RPolyCtrl\n";
  
  numSubjectsUsed = geneIndex.size();
  numCov = Rcpp::as<std::vector<double> >(subjectData["covariates"]).size() / numSubjectsUsed;
  ConvertCovariates(subjectData, covariates, numCov, numSubjectsUsed);
//  Rcpp::Rcout << "NUmber of subjects used:\t" << numSubjectsUsed << std::endl;
//  Rcpp::Rcout << "Number of covariate values:\t" << numCov << std::endl;
  
  geneticData = OpenGeneticData(geneticInfo);
  if (geneticData == NULL)
    Rcpp::stop("Unable to open genetic data");
  if (geneticData->GetFirst()) {
    Rcpp::Rcout << "Error reading first SNP" << std::endl;
    delete geneticData;
    outfile.close();
    return 1;
  }
//  geneticData->GetNext();
//  geneticData->GetNext();
  
  // Lot of legacy code. There should be no missing covariates or phenotypes
  missingPhenotype.assign(phenotype.size(), false);
  caseControl.assign(phenotype.size(), false);
  for (ui = 0; ui < phenotype.size(); ++ui) {
    if (phenotype[ui] == 1)
      caseControl[ui] = true;
    else if (phenotype[ui] != 0)
      missingPhenotype[ui] = true;
  }
  missingCov.assign(covariates.size(), false);
  missingGene.assign(numSubjectsUsed, false); // This should be returned from CGeneticData
  filter.assign(numSubjectsUsed, true); // All are true - no filtering is currently supported
  if (geneticData->GeneticProbabilities())
    geneticValues.resize(4 * numSubjectsUsed);
  else
    geneticValues.resize(numSubjectsUsed);
  CGxEPolytomousDataset gxeData(numSubjectsUsed, caseControl.data(), (bool *)missingPhenotype.data(),
                                numCov, covariates.data(), (bool *)missingCov.data(), numCov, (bool *)filter.data());

  if (gxeData.Initialize() == false)
    Rcpp::Rcout << "Initialization failed - " << gxeData.ErrorString() << std::endl;
  
  d = &geneticValues[0];
  if (geneticData->GeneticProbabilities()) {
    p0 = d + numSubjectsUsed;
    p1 = p0 + numSubjectsUsed;
    p2 = p1 + numSubjectsUsed;
    for (i = 0; i < numSubjectsUsed; ++i, ++d, ++p0, ++p1, ++p2) {
//      if (i < 5)
//        Rcpp::Rcout << geneIndex[i] << '\t';
      *d = geneticData->Dosages()[geneIndex[i] - 1];
      *p0 = geneticData->Probabilities()[0][geneIndex[i] - 1];
      *p1 = geneticData->Probabilities()[1][geneIndex[i] - 1];
      *p2 = geneticData->Probabilities()[2][geneIndex[i] - 1];
    }
//    Rcpp::Rcout << std::endl;
  } else {
    for (i = 0; i < numSubjectsUsed; ++i, ++d)
      *d = geneticData->Dosages()[geneIndex[i] - 1];
  }
  
  if (geneticData->GeneticProbabilities())
    gxeData.AssignGene(geneticValues.data(), (bool *)missingGene.data(), true, true);
  else
    gxeData.AssignGene(geneticValues.data(), (bool *)missingGene.data(), true, false);
  gxeData.UpdateGene();
//  gxeData.Print(outfile);
  
  retVal = gxeData.FitModels();
  if (chrName[0] != "")
    outfile << chrName[0] << '\t';
  if (snpName[0] != "")
    outfile << snpName[0] << '\t';
  if (bp[0] != 0)
    outfile << bp[0] << '\t';
  if (a1Name[0] != "") {
    if (gxeData.AllelesSwapped())
      outfile << a2Name[0] << '\t' << a1Name[0] << '\t';
    else
      outfile << a1Name[0] << '\t' << a2Name[0] << '\t';
  }
  outfile << gxeData.NumCasesUsed() << '\t' << gxeData.NumControlsUsed() << '\t';
  WriteResults(outfile, retVal, gxeData);

  for (j = 1; j < geneticData->NumSNPs(); ++j) {
    geneticData->GetNext();
    d = &geneticValues[0];
    if (geneticData->GeneticProbabilities()) {
      p0 = d + numSubjectsUsed;
      p1 = p0 + numSubjectsUsed;
      p2 = p1 + numSubjectsUsed;
      for (i = 0; i < numSubjectsUsed; ++i, ++d, ++p0, ++p1, ++p2) {
//        if (i < 5)
//          Rcpp::Rcout << geneIndex[i] << '\t';
        *d = geneticData->Dosages()[geneIndex[i] - 1];
        *p0 = geneticData->Probabilities()[0][geneIndex[i] - 1];
        *p1 = geneticData->Probabilities()[1][geneIndex[i] - 1];
        *p2 = geneticData->Probabilities()[2][geneIndex[i] - 1];
      }
//      Rcpp::Rcout << std::endl;
    } else {
      for (i = 0; i < numSubjectsUsed; ++i, ++d)
        *d = geneticData->Dosages()[geneIndex[i] - 1];
    }
    gxeData.UpdateGene();
    retVal = gxeData.FitModels();
    if (chrName[0] != "")
      outfile << chrName[j] << '\t';
    if (snpName[0] != "")
      outfile << snpName[j] << '\t';
    if (bp[0] != 0)
      outfile << bp[j] << '\t';
    if (a1Name[0] != "") {
      if (gxeData.AllelesSwapped())
        outfile << a2Name[j] << '\t' << a1Name[j] << '\t';
      else
        outfile << a1Name[j] << '\t' << a2Name[j] << '\t';
    }
    outfile << gxeData.NumCasesUsed() << '\t' << gxeData.NumControlsUsed() << '\t';
    WriteResults(outfile, retVal, gxeData);
  }  
  
  if (geneticData)
    delete geneticData;
  
/*    
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
  j = 1;
  while(geneticData->GetNext() == 0 && j < 5) {
    ++j;
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
    gxeData.UpdateGene();
    retVal = gxeData.FitModels();
//    Rcpp::Rcout << "Fit Model, result:\t" << retVal << std::endl;
    WriteResults(outfile, retVal, gxeData);
  }
  
  if (geneticData)
    delete geneticData;

  std::vector<double> covInt;
  covInt.assign(numSubjectsUsed * gxeData.NumParameters(), -9);
  if (gxeData.CompleteCovariates() != NULL) {
    for (i = 0; i < numSubjectsUsed * gxeData.NumParameters(); ++i)
      covInt[i] = gxeData.CompleteCovariates()[i];
  }
  
  std::ofstream outfile("GxEDataset.txt");
  for (i = 0; i < numSubjects; ++i) {
    outfile << phenotype[i];
    for (j = 0; j < gxeData.NumParameters(); ++j)
      outfile << '\t' << covInt[i * gxeData.NumParameters() + j];
    outfile << std::endl;
  }
  outfile.close();
 
//  Rcpp::Rcout << "Number of parameters\t" << gxeData.NumParameters() << std::endl;
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
 */
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