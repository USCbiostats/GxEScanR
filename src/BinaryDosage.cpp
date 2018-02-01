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
Rcpp::List BinaryDosageInfo(const std::string &geneticFilename, const std::string &mapFilename, const std::string &familyFilename) {
  Rcpp::List result;
  std::string errorMessage = "";
  // Needed files, genetic file (always needed), map and family file (optional)
  Rcpp::StringVector files(3);
  // Version and subversion number
  Rcpp::IntegerVector version(2);
  // Number of subjects, SNPs, and groups
  Rcpp::IntegerVector counts(3);
  unsigned int numSubjects, numSNPs, numGroups;
  
  files[0] = geneticFilename;
  files[1] = mapFilename;
  files[2] = familyFilename;
  
  result = Rcpp::List::create(Rcpp::Named("type") = "binaryDosage",
                              Rcpp::Named("files") = files,
                              Rcpp::Named("valid") = false,
                              Rcpp::Named("errorMessage") = "",
                              Rcpp::Named("Version") = version,
                              Rcpp::Named("Counts") = counts);
  // Make sure the genetic file is a binary dosage file and get the version number
  if (GetBinaryDosageFormat(geneticFilename, version[0], version[1], errorMessage)) {
    result["errorMessage"] = errorMessage;
    return result;
  }
  // Get the number of Subjects, SNPs, and groups
  if (version[0] < 4) {
    // For versions 1, 2, and 3, the number os subjects comes from the family file, and
    // the number of SNPs comes from the map file, and the number of groups is 1.
    if (GetNumberOfSubjectsAndSNPs(familyFilename, mapFilename, numSubjects, numSNPs, errorMessage)) {
      result["errorMessage"] = errorMessage;
      return result;
    }
    counts[0] = numSubjects;
    counts[1] = numSNPs;
    counts[2] = 1;
    // For versions 1 and 2, the file size is fixed based on the number of subjects and SNPs.
    if (version[0] < 3) {
      if (CheckBinaryDosageFileSize(geneticFilename, version[1], numSubjects, numSNPs, errorMessage)) {
        result["errorMessage"] = errorMessage;
        return result;
      }
    }
  } else {
    // For version 4 the number of subjects, SNPs, and groups are in the binary dosage file
    if (GetNumberOfSubjectsAndSNPsV4(geneticFilename, numSubjects, numSNPs, numGroups, errorMessage)) {
      result["errorMessage"] = errorMessage;
      return result;
    }
    counts[0] = numSubjects;
    counts[1] = numSNPs;
    counts[2] = numGroups;
  }

  result["valid"] = true;
  return result;
}

int GetBinaryDosageFormat(const std::string &binaryDosageFilename, int &format, int &version, std::string &errorMessage) {
  std::ifstream infile;
  const char magicWord[4] = { 'b', 'o', 's', 'e' };
  char readMagicWord[4];
  char readVersion[4];
  
  format = 0;
  version = 0;
  
  infile.open(binaryDosageFilename.c_str(), std::ios_base::in | std::ios_base::binary);
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

int GetNumberOfSubjectsAndSNPs(const std::string &familyFilename, const std::string &mapFilename, unsigned int &numSubjects, unsigned int &numSNPs, std::string &errorMessage) {
  Rcpp::Environment env = Rcpp::Environment::namespace_env("GxEScanR");
  Rcpp::Function testFamFile = env["FamilyFileCheck"];
  Rcpp::Function testMapFile = env["MapFileCheck"];
  Rcpp::NumericVector checkResult = testFamFile(familyFilename);
  
  if (checkResult[0] == 0) {
    errorMessage = "Error reading family file";
    return 1;
  }
  numSubjects = checkResult[0];
  
  checkResult = testMapFile(mapFilename);
  if (checkResult[0] == 0) {
    errorMessage = "Error reading map file";
    return 1;
  }
  numSNPs = checkResult[0];
  return 0;
}

int GetNumberOfSubjectsAndSNPsV4(const std::string &geneticFilename, unsigned int &numSubjects, unsigned int &numSNPs, unsigned int &numGroups, std::string &errorMessage) {
  std::ifstream infile;
  int counts[3];
  
  infile.open(geneticFilename);
  if (!infile.good()) {
    infile.close();
    errorMessage = "Unable to open binary dosage file";
    return 1;
  }
  infile.seekg(8);
  infile.read((char *)&counts[0], 3 * sizeof(unsigned int));
  if (!infile.good()) {
    infile.close();
    errorMessage = "Error reading number of subjects, SNPs, and groups from binary dosage file";
    return 1;
  }
  numSubjects = counts[0];
  numSNPs = counts[1];
  numGroups = counts[2];
  if (numSubjects == 0) {
    errorMessage = "Number of subjects in binary dosage file equals 0";
    return 1;
  }
  if (numSNPs == 0) {
    errorMessage = "Number of SNPs in binary dosage file equals 0";
    return 1;
  }
  if (numGroups == 0) {
    errorMessage = "Number of groups in binary dosage file equals 0";
    return 1;
  }
  
  infile.close();
  return 0;
}

// Only applies to formats 1 and 2
int CheckBinaryDosageFileSize(const std::string &geneticFilename, int version, unsigned int numSubjects, unsigned int numSNPs, std::string &errorMessage) {
  // The expected file size
  long long expectedFileSize;
  // Actual size of file
  std::streampos actualFileSize;
  std::ifstream infile;
  
  expectedFileSize = sizeof(unsigned short) * numSubjects * numSNPs * version + 8;
  infile.open(geneticFilename);
  infile.seekg(0, std::ios_base::end);
  actualFileSize = infile.tellg();
  infile.close();
  if (actualFileSize != expectedFileSize) {
    errorMessage = "Binary dosage file is not the correct size";
    return 1;
  }
  return 0;
}

CBinaryDosage::CBinaryDosage() : CGeneticData() {
  m_geneticFilename = "";
  m_version = 0;
  m_subversion = 0;
}

CBinaryDosage::CBinaryDosage(std::string &_geneticFilename) : CGeneticData() {
  m_geneticFilename = _geneticFilename;
  m_version = 0;
  m_subversion = 0;
}

CBinaryDosageFormat1::CBinaryDosageFormat1(std::string &_geneticFilename, std::string &_mapFilename, std::string &_familyFilename) : CBinaryDosage(_geneticFilename) {
  m_mapFilename = _mapFilename;
  m_familyFilename = _familyFilename;
}

CBinaryDosageFormat1::CBinaryDosageFormat1(Rcpp::List &_binaryDosageInfo) : CBinaryDosage() {

  Rcpp::StringVector fileVector;
  Rcpp::IntegerVector versionVector;
  Rcpp::IntegerVector countVector;

  fileVector = Rcpp::as<Rcpp::StringVector>(_binaryDosageInfo["files"]);
  m_geneticFilename = Rcpp::as<std::string>(fileVector[0]);
  m_mapFilename = Rcpp::as<std::string>(fileVector[1]);
  m_familyFilename = Rcpp::as<std::string>(fileVector[2]);
 
  versionVector = Rcpp::as<Rcpp::IntegerVector>(_binaryDosageInfo["Version"]);
  m_version = versionVector[0];
  m_subversion = versionVector[1];
  
  countVector = Rcpp::as<Rcpp::IntegerVector>(_binaryDosageInfo["Counts"]);
  m_numSubjects = countVector[0];
  m_numSNPs = countVector[1];
}

int CBinaryDosageFormat1::CheckValidity() {
  if (m_version == 0 && m_subversion == 0) {
    if (GetBinaryDosageFormat(m_geneticFilename, m_version, m_subversion, m_errorMessage))
      return 1;
  }
  if (CheckVersion())
    return 1;
  
  if (m_numSubjects != 0 && m_numSNPs != 0) {
    if (GetNumberOfSubjectsAndSNPs(m_familyFilename, m_mapFilename, m_numSubjects, m_numSNPs, m_errorMessage))
      return 1;
  }
  if (CheckBinaryDosageFileSize(m_geneticFilename, m_subversion, m_numSubjects, m_numSNPs, m_errorMessage))
    return 1;
  
  m_valid = true;
  return 0;
}

CBinaryDosageFormat1_1::CBinaryDosageFormat1_1(std::string &_geneticFilename, std::string &_mapFilename, std::string &_familyFilename) :
  CBinaryDosageFormat1(_geneticFilename, _mapFilename, _familyFilename) {
}

CBinaryDosageFormat1_1::CBinaryDosageFormat1_1(Rcpp::List &_binaryDosageInfo) : CBinaryDosageFormat1(_binaryDosageInfo) {
}

int CBinaryDosageFormat1_1::CheckVersion() {
  if (m_version != 1 || m_subversion != 1) {
    m_errorMessage = "Binary dosage is not in format 1.1";
    return 1;
  }
  return 0;
}
// gfn <- "C:/GxEScan/Hita/bchr21.bdosage"
// mfn <- "C:/GxEScan/Hita/bchr21.bim"
// ffn <- "C:/GxEScan/Hita/chs3000.fam"
// pbi <- BinaryDosageInfo(gfn, mfn, ffn)
// GxEScan(df2, pbi)
