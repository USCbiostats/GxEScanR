#include <string>
#include <RcppArmadillo.h>
#include "GeneticData.h"
#include "PlinkBinary.h"

const unsigned int NumGeneticDataTypes = 1;
const std::string GeneticDataTypes[NumGeneticDataTypes] = { "plinkBinary" };

CGeneticData *TestGeneticData(Rcpp::List &_geneticData) {
  CGeneticData *geneticData = NULL;
  Rcpp::StringVector typeVector;
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
    if (geneticData->CheckValidity()) {
      delete geneticData;
      geneticData = NULL;
    }
  default:
    break;
  }
  return geneticData;
}

CGeneticData::CGeneticData() {
  m_errorMessage = "";
  m_numSubjects = 0;
  m_numSNPs = 0;
}
