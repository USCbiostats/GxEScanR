#include <string>
#include <RcppArmadillo.h>
#include "GeneticData.h"
#include "PlinkBinary.h"
#include "BinaryDosage.h"

const unsigned int NumGeneticDataTypes = 2;
const std::string GeneticDataTypes[NumGeneticDataTypes] = { "plinkBinary", "binaryDosage" };

CGeneticData *TestGeneticData(Rcpp::List &_geneticData) {
  CGeneticData *geneticData = NULL;
  Rcpp::StringVector typeVector;
  Rcpp::IntegerVector versionVector;
  std::string type;
  unsigned int ui;
  
  typeVector = Rcpp::as<Rcpp::StringVector>(_geneticData["type"]);
  type = Rcpp::as<std::string>(typeVector[0]);
  for (ui = 0; ui < NumGeneticDataTypes; ++ui) {
    if (type == GeneticDataTypes[ui])
      break;
  }
  switch(ui) {
  case 0:
    geneticData = new CPlinkBinary(_geneticData);
    break;
  case 1:
    versionVector = Rcpp::as<Rcpp::IntegerVector>(_geneticData["Version"]);
    switch(versionVector[0]) {
    case 1:
      switch(versionVector[1]) {
      case 1:
        Rcpp::Rcout << "Before creating BinaryDosageFormat1_1" << std::endl;
        geneticData = new CBinaryDosageFormat1_1(_geneticData);
        break;
      default:
        geneticData = NULL;
        break;
      }
      break;
    default:
      geneticData = NULL;
      break;
    }
    break;
  default:
    geneticData = NULL;
    break;
  }
  if (geneticData == NULL)
    return geneticData;
  if (geneticData->CheckValidity()) {
    delete geneticData;
    geneticData = NULL;
  }
  return geneticData;
}

CGeneticData::CGeneticData() {
  m_errorMessage = "";
  m_numSubjects = 0;
  m_numSNPs = 0;
  m_valid = false;
  m_bMeasured = false;
  m_bGeneticProbabilities = false;
}

int CGeneticData::GetSNP(unsigned int n) {
  if (!m_valid) {
    m_errorMessage = "Data for genetic file is invalid";
    return 1;
  }
  if (n == 0 || n > m_numSNPs) {
    m_errorMessage = "Invalid SNP value";
    return 1;
  }
  return 0;
}