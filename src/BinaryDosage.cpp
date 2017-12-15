#include <istream>
#include <fstream>
#include <string>
#include <RcppArmadillo.h>
#include "BinaryDosage.h"

//' Function to produce summary of a binary dosage file used by GxEScan
//' 
//' Function to produce summary of a binary dosage file used by GxEScan
//' 
//' @param geneticFile
//' Name of file with genetic data, normally ends with .bdosage
//' @param mapFile
//' Name of map file associated with genetic data file. Not needed for format version 4.
//' The map is the exteneded plink format and normally ends with .bim
//' @param familyFile
//' Name of family file associated with genetic data file. Not needed for format version 4.
//' The family file is the plink family fiel format and normally ends with .fam
//' @return
//' List with data needed by GxEScan
//' @export
// [[Rcpp::export]]
Rcpp::List BinaryDosageInfo(std::string &geneticFile, std::string &mapFile, std::string &familyFile) {
  Rcpp::List result;
  std::string errorMessage = "";
  // Needed files, genetic file (always needed), map and family file (optional)
  Rcpp::StringVector files(3);
  // Version and subversion number
  Rcpp::IntegerVector version(2);
  // Number of subjects, SNPs, and groups
  Rcpp::IntegerVector counts(3);
  // The expected file size - only needed for formats 1.x and 2.x.
  long long expectedFileSize;
  // Actual size of file
  std::streampos actualFileSize;
  std::ifstream infile;
  
  files[0] = geneticFile;
  files[1] = mapFile;
  files[2] = familyFile;
  
  result = Rcpp::List::create(Rcpp::Named("type") = "binaryDosage",
                              Rcpp::Named("files") = files,
                              Rcpp::Named("valid") = false,
                              Rcpp::Named("errorMessage") = "",
                              Rcpp::Named("Version") = version,
                              Rcpp::Named("Counts") = counts);

  if (GetBinaryDosageFormat(geneticFile, version[0], version[1], errorMessage)) {
    result["errorMessage"] = errorMessage;
    return result;
  }

  if (version[0] < 4) {
    Rcpp::Environment env = Rcpp::Environment::namespace_env("GxEScanR");
    Rcpp::Function testFamFile = env["FamilyFileCheck"];
    Rcpp::Function testMapFile = env["MapFileCheck"];
    Rcpp::NumericVector checkResult = testFamFile(familyFile);
    
    if (checkResult[0] == 0) {
      result["errorMessage"] = "Error reading family file";
      return result;
    }
    counts[0] = checkResult[0];

    checkResult = testMapFile(mapFile);
    if (checkResult[0] == 0) {
      result["errorMessage"] = "Error reading map file";
      return result;
    }
    counts[1] = checkResult[0];
  } else {
    infile.open(geneticFile);
    infile.seekg(8);
    infile.read((char *)&counts[0], 3 * sizeof(unsigned int));
    if (!infile.good()) {
      infile.close();
      result["errorMessage"] = "Error reading number of subjects, SNPs, and groups from binary dosage file";
      return result;
    }
    infile.close();
  }
  
  if (version[0] < 3) {
    expectedFileSize = 2 * counts[0] * counts[1] * version[1] + 8;
    infile.open(geneticFile);
    infile.seekg(0, std::ios_base::end);
    actualFileSize = infile.tellg();
    infile.close();
    if (actualFileSize != expectedFileSize) {
      result["errorMessage"] = "Binary dosage file is not the correct size";
      return result;
    }
    counts[2] = 1;
  }

  result["valid"] = true;
  return result;
}

int GetBinaryDosageFormat(std::string &binaryDosageFile, int &format, int &version, std::string &errorMessage) {
  std::ifstream infile;
  const char magicWord[4] = { 'b', 'o', 's', 'e' };
  char readMagicWord[4];
  char readVersion[4];
  
  format = 0;
  version = 0;
  
  infile.open(binaryDosageFile.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!infile.good()) {
    errorMessage = "Unable to open binary dosage file";
    return 1;
  }
  infile.read(readMagicWord, 4);
  if (!infile.good()) {
    infile.close();
    errorMessage = "Error reading header of binary dosage file";
    return 1;
  }
  if (std::memcmp(magicWord, readMagicWord, 4)) {
    infile.close();
    errorMessage = "File is not a binary dosage file";
    return 1;
  }
  infile.read(readVersion, 4);
  if (!infile.good()) {
    infile.close();
    errorMessage = "Error reading format and version number";
    return 1;
  }
  infile.close();
  if (readVersion[0] || readVersion[2]) {
    errorMessage = "Format and version number error";
    return 1;
  }
  format = readVersion[1];
  version = readVersion[3];
  if (format < 1 || format > 4) {
    errorMessage = "Unknown format";
    return 1;
  }
  if (version < 1 || version > 2) {
    errorMessage = "Unknown version";
    return 1;
  }
  
  infile.close();
  return 0;
}

// gfn <- "C:/GxEScan/Hita/bchr21.bdosage"
// mfn <- "C:/GxEScan/Hita/bchr21.bim"
// ffn <- "C:/GxEScan/Hita/chs3000.fam"
// pbi <- BinaryDosageInfo(gfn, mfn, ffn)
// GxEScan(df2, pbi)

CBinaryDosage::CBinaryDosage() : CGeneticData() {
  m_geneticFile = "";
  m_version = 0;
  m_subversion = 0;
}

CBinaryDosage::CBinaryDosage(std::string &_geneticFile) : CGeneticData() {
  m_geneticFile = _geneticFile;
  m_version = 0;
  m_subversion = 0;
}

CBinaryDosageFormat1::CBinaryDosageFormat1(std::string &_geneticFile, std::string &_mapFile, std::string &_familyFile) : CBinaryDosage(_geneticFile) {
  m_mapFile = _mapFile;
  m_familyFile = _familyFile;
}

CBinaryDosageFormat1::CBinaryDosageFormat1(Rcpp::List &_binaryDosageInfo) : CBinaryDosage() {

  Rcpp::StringVector fileVector;
  Rcpp::IntegerVector versionVector;
  Rcpp::IntegerVector countVector;

  fileVector = Rcpp::as<Rcpp::StringVector>(_binaryDosageInfo["files"]);
  m_geneticFile = Rcpp::as<std::string>(fileVector[0]);
  m_mapFile = Rcpp::as<std::string>(fileVector[1]);
  m_familyFile = Rcpp::as<std::string>(fileVector[2]);
 
  versionVector = Rcpp::as<Rcpp::IntegerVector>(_binaryDosageInfo["Version"]);
  m_version = versionVector[0];
  m_subversion = versionVector[1];
  
  countVector = Rcpp::as<Rcpp::IntegerVector>(_binaryDosageInfo["Counts"]);
  m_numSubjects = countVector[0];
  m_numSNPs = countVector[1];
}

CBinaryDosageFormat1_1::CBinaryDosageFormat1_1(std::string &_geneticFile, std::string &_mapFile, std::string &_familyFile) :
  CBinaryDosageFormat1(_geneticFile, _mapFile, _familyFile) {
}

CBinaryDosageFormat1_1::CBinaryDosageFormat1_1(Rcpp::List &_binaryDosageInfo) : CBinaryDosageFormat1(_binaryDosageInfo) {
  
}