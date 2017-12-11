#include <istream>
#include <fstream>
#include <string>
#include <RcppArmadillo.h>
#include "BinaryDosage.h"

//' Function to produce summary of a plink binary genetic file used by GxEScan
//' 
//' Function to produce summary of a plink binary genetic file used by GxEScan
//' 
//' @param geneticFile
//' Name of file with genetic data, normally ends with .bdosage
//' @param mapFile
//' Name of map file associated with genetic data file.
//' The map is the exteneded plink format and normally ends with .bim
//' @param familyFile
//' Name of family file associated with genetic data file.
//' The family file is the plink family fiel format and normally ends with .fam
//' @return
//' List with data needed by GxEScan
//' @export
// [[Rcpp::export]]
Rcpp::List BinaryDosageInfo(std::string &geneticFile, std::string &mapFile, std::string &familyFile) {
  Rcpp::List result;
  std::string errorMessage = "";
  Rcpp::StringVector files(3);
  int format, version;
  int numSubjects, numSNPs;
  long long expectedFileSize;
  std::streampos actualFileSize;
  std::ifstream infile;
  
  files[0] = geneticFile;
  files[1] = mapFile;
  files[2] = familyFile;
  
  result = Rcpp::List::create(Rcpp::Named("type") = "binaryDosage",
                              Rcpp::Named("files") = files,
                              Rcpp::Named("valid") = false,
                              Rcpp::Named("errorMessage") = errorMessage,
                              Rcpp::Named("Format") = 0,
                              Rcpp::Named("Version") = 0,
                              Rcpp::Named("numSubjects") = 0,
                              Rcpp::Named("numSNPs") = 0);

  if (GetBinaryDosageFormat(geneticFile, format, version, errorMessage)) {
    result["errorMessage"] = errorMessage;
    return result;
  }
  
  result["Format"] = format;
  result["Version"] = version;
  if (format < 4) {
    Rcpp::Environment env = Rcpp::Environment::namespace_env("GxEScanR");
    Rcpp::Function testFamFile = env["FamilyFileCheck"];
    Rcpp::Function testMapFile = env["MapFileCheck"];
    Rcpp::NumericVector checkResult = testFamFile(familyFile);
    
    if (checkResult[0] == 0) {
      errorMessage = "Error reading family file";
      result["errorMessage"] = errorMessage;
      return result;
    }
    numSubjects = checkResult[0];
    result["numSubjects"] = checkResult;

    checkResult = testMapFile(mapFile);
    if (checkResult[0] == 0) {
      errorMessage = "Error reading map file";
      result["errorMessage"] = errorMessage;
      return result;
    }
    numSNPs = checkResult[0];
    result["numSNPs"] = checkResult;
  }
  
  if (format < 3) {
    expectedFileSize = 2 * numSubjects * numSNPs * version + 8;
    infile.open(geneticFile);
    infile.seekg(0, std::ios_base::end);
    actualFileSize = infile.tellg();
    infile.close();
    if (actualFileSize != expectedFileSize) {
      errorMessage = "Binary dosage file is not the correct size";
      return result;
    }
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