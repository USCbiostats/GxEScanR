#include <string>
#include <RcppArmadillo.h>
#include "GeneticData.h"

const unsigned int NumGeneticDataTypes = 1;
const std::string GeneticDataTypes[NumGeneticDataTypes] = { "plinkBinary" };

int TestGeneticData(Rcpp::List &_geneticData) {
  unsigned int ui;
  Rcpp::StringVector typeVector;
  std::string type;
  
  typeVector = Rcpp::as<Rcpp::StringVector>(_geneticData["type"]);
  type = Rcpp::as<std::string>(typeVector[0]);
  for (ui = 0; ui < NumGeneticDataTypes; ++ui) {
    if (type == GeneticDataTypes[ui])
      break;
  }
  if (ui == NumGeneticDataTypes) {
    Rcpp::Rcout << "Unknown genetic data type" << std::endl;
    return 1;
  }
  return 0;
}

CGeneticData::CGeneticData() {
  m_errorMessage = "";
  m_numSubjects = 0;
  m_numSNPs = 0;
}
