#include <istream>
#include <fstream>
#include <RcppArmadillo.h>
#include "PlinkBinary.h"

//' Function to produce summary of a plink binary genetic file used by GxEScan
//' 
//' Function to produce summary of a plink binary genetic file used by GxEScan
//' 
//' @param geneticFile
//' Name of file with genetic data, normally ends with .bed
//' @param mapFile
//' Name of map file associated with genetic data file, normally ends with .bim
//' @param familyFile
//' Name of family file associated with genetic data file, normally ends with .fam
//' @return
//' List with data needed by GxEScan
//' @export
// [[Rcpp::export]]
Rcpp::List PlinkBinaryInfo(std::string &geneticFile, std::string &mapFile, std::string &familyFile) {
  Rcpp::List result;
  std::string type = "plinkBinary";
  CPlinkBinary geneticData(geneticFile, mapFile, familyFile);
  
  result = Rcpp::List::create(Rcpp::Named("type") = type,
                              Rcpp::Named("file") = geneticFile,
                              Rcpp::Named("mapFile") = mapFile,
                              Rcpp::Named("familyFile") = familyFile);
  
  if (geneticData.CheckValidity()) {
    Rcpp::Rcout << geneticData.ErrorMessage() << std::endl;
  }
  
  return result;
}

CPlinkBinary::CPlinkBinary(std::string &_geneticFile, std::string &_mapFile, std::string &_familyFile) : CGeneticData() {
  m_geneticFile = _geneticFile;
  m_mapFile = _mapFile;
  m_familyFile = _familyFile;
}

int CPlinkBinary::CheckValidity() {
  std::ifstream infile;
  
  infile.open(m_geneticFile.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!infile.good()) {
    m_errorMessage = "Unable to open plink binary genetic file";
    return 1;
  }
  infile.close();
  infile.open(m_mapFile.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!infile.good()) {
    m_errorMessage = "Unable to open map file";
    return 1;
  }
  infile.close();
  infile.open(m_familyFile.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!infile.good()) {
    m_errorMessage = "Unable to open family file";
    return 1;
  }
  infile.close();
  
  Rcpp::Environment env = Rcpp::Environment::namespace_env("GxEScanR");
  Rcpp::Function testSubjectFunction = env["FamilyFileInfo"];
  Rcpp::NumericVector result = testSubjectFunction(m_familyFile);
  
  if (result[0] == 0) {
    Rcpp::Rcout << "Error reading family file" << std::endl;
    return 1;
  }
  Rcpp::Rcout << result[0] << " subjects in family file" << std::endl;
  return 0;
}

// gfn <- "C:/GxEScan/SmallData.bed"
// mfn <- "C:/GxEScan/SmallData.bim"
// ffn <- "C:/GxEScan/SmallData.fam"
// pbi <- PlinkBinaryInfo(gfn, mfn, ffn)