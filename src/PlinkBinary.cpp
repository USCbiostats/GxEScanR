#include <istream>
#include <fstream>
#include <string>
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
  CPlinkBinary geneticData(geneticFile, mapFile, familyFile);
  Rcpp::StringVector files(3);
  
  files[0] = geneticFile;
  files[1] = mapFile;
  files[2] = familyFile;
  
  result = Rcpp::List::create(Rcpp::Named("type") = "plinkBinary",
                              Rcpp::Named("files") = files,
                              Rcpp::Named("valid") = false,
                              Rcpp::Named("errorMessage") = "",
                              Rcpp::Named("numSubjects") = 0,
                              Rcpp::Named("numSNPs") = 0);
  
  if (geneticData.CheckValidity()) {
    result["errorMessage"] = geneticData.ErrorMessage();
  } else {
    result["valid"] = true;
    result["numSubjects"] = geneticData.NumSubjects();
    result["numSNPs"] = geneticData.NumSNPs();
  }
  
  return result;
}

CPlinkBinary::CPlinkBinary(Rcpp::List &_geneticData) {
  Rcpp::StringVector fileVector;

  fileVector = Rcpp::as<Rcpp::StringVector>(_geneticData["files"]);
  m_geneticFile = Rcpp::as<std::string>(fileVector[0]);
  m_mapFile = Rcpp::as<std::string>(fileVector[1]);
  m_familyFile = Rcpp::as<std::string>(fileVector[2]);
}

CPlinkBinary::CPlinkBinary(std::string &_geneticFile, std::string &_mapFile, std::string &_familyFile) : CGeneticData() {
  m_geneticFile = _geneticFile;
  m_mapFile = _mapFile;
  m_familyFile = _familyFile;
}

int CPlinkBinary::CheckFileSize() {
  long long expectedFileSize;
  std::streampos actualFileSize;
  std::ifstream infile;
  char expectedHeader[3] = { 0x6c, 0x1b, 0x01 };
  char header[3];
  
  infile.open(m_geneticFile, std::ios_base::in | std::ios_base::binary);
  infile.read(header, 3);
  if (!infile.good()) {
    Rcpp::Rcout << "Error reading header" << std::endl;
    infile.close();
    return 1;
  }
  if (std::memcmp(header, expectedHeader, 3)) {
    Rcpp::Rcout << "Genetic file is not a plink binary formatted file" << std::endl;
    infile.close();
    return 1;
  }
  infile.seekg(0, std::ios_base::end);
  actualFileSize = infile.tellg();
  infile.close();
  
  expectedFileSize = ((m_numSubjects + 3) / 4) * m_numSNPs + 3;
  if (expectedFileSize != actualFileSize) {
    Rcpp::Rcout << "Expected file size: " << expectedFileSize << ", does not equal actual file size: " << actualFileSize << std::endl;
    return 1;
  }
  return 0;
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
  Rcpp::Function testFamFile = env["FamilyFileCheck"];
  Rcpp::NumericVector result = testFamFile(m_familyFile);

  if (result[0] == 0) {
    Rcpp::Rcout << "Error reading family file" << std::endl;
    return 1;
  }
  m_numSubjects = result[0];

  Rcpp::Function testMapFile = env["MapFileCheck"];
  result = testMapFile(m_mapFile);
  
  if (result[0] == 0) {
    Rcpp::Rcout << "Error reading map file" << std::endl;
    return 1;
  }
  m_numSNPs = result[0];

  if (CheckFileSize())
    return 1;
  
  m_valid = true;
  return 0;
}

// gfn <- "C:/GxEScan/SmallData.bed"
// mfn <- "C:/GxEScan/SmallData.bim"
// ffn <- "C:/GxEScan/SmallData.fam"
// pbi <- PlinkBinaryInfo(gfn, mfn, ffn)
// GxEScan(df2, pbi)