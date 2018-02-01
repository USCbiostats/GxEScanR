#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <Rcpp.h>
using namespace Rcpp;

bool CheckBinaryDosageFileSize(std::ifstream &infile, const int version, const unsigned int nSub, const unsigned int nSNPs) {
  long long expectedSize;
  std::streampos actualSize;
  
  expectedSize = nSub * nSNPs * sizeof(short);
  if (version == 2)
    expectedSize += expectedSize;
  expectedSize += 8;
  infile.seekg(0, std::ios_base::end);
  actualSize = infile.tellg();
  Rcpp::Rcout << "Format:\tX." << version << '\t' << expectedSize << '\t' << actualSize << std::endl;
  if (expectedSize == actualSize)
    return true;
  Rcpp::Rcerr << "Binary dosage file is not of the expected size" << std::endl;
  return false;
}

bool CheckBinaryDosage3FileSize(std::ifstream &infile, const int version, const unsigned int nSub, const unsigned int nSNPs) {
  long long expectedSize;
  int snpSize;
  std::streampos actualSize;
  unsigned int ui;
  
  
  if (version == 1) {
    expectedSize = nSub * nSNPs * sizeof(short) + 12;
    infile.seekg(0, std::ios_base::end);
    actualSize = infile.tellg();
    Rcpp::Rcout << "Format:\t3.1\t" << expectedSize << '\t' << actualSize << std::endl;
    if (expectedSize == actualSize)
      return true;
    Rcpp::Rcerr << "Binary dosage file is not of the expected size" << std::endl;
    return false;
  } else {
    infile.seekg(0, std::ios_base::end);
    actualSize = infile.tellg();
    expectedSize = 12 + 4 * nSNPs;
    infile.seekg(12);
    for (ui = 0; ui < nSNPs; ++ui) {
      infile.read((char *)&snpSize, sizeof(int));
      if (!infile.good()) {
        Rcpp::Rcerr << "Reached end of file:\t" << expectedSize << '\t' << actualSize << '\t' << snpSize << '\t' << infile.tellg() << std::endl;
        return false;
      }
      if (ui < 10)
        Rcpp::Rcout << snpSize << std::endl;
      expectedSize += snpSize;
      infile.seekg(snpSize, std::ios_base::cur);
    }
    if (expectedSize != actualSize) {
      Rcpp::Rcerr << "Binary dosage file is not of the expected size\t" << expectedSize << '\t' << actualSize << std::endl;
      return false;
    }
    Rcpp::Rcout << "Format:\t3.1\t" << expectedSize << '\t' << actualSize << std::endl;
  }
  
  return true;
}

void ReadStringArray(std::ifstream &infile, std::vector<std::string> &x, unsigned int numStrings, unsigned int arraySize) {
  char *charArray = NULL;
  std::string stringArray;
  unsigned int ui;

  x.resize(numStrings);
  charArray = new char[arraySize];
  infile.read(charArray, arraySize);
  stringArray = charArray;

  std::istringstream is(stringArray);
  for (ui = 0; ui < numStrings; ++ui)
    is >> x[ui];
  
  if (charArray)
    delete [] charArray;
  
}

bool GetBinaryDosage4Info(std::ifstream &infile, unsigned int &numSubjects, unsigned int &numSNPs, Rcpp::List &result) {
  Rcpp::List subjects;
  Rcpp::List snps;
  int numGroups;
  unsigned int subjectOptions;
  unsigned int snpOptions;
  unsigned int uiStartSubjects, uiStartSNPs, uiStartDosage;
  std::streampos startSubjects, startSNPs, startDosage;
  unsigned int subjectStringSize, familyStringSize, snpStringSize, chromosomeStringSize, refAlleleStringSize, altAlleleStringSize;
  std::vector<std::string> subjectID, familyID, snpID, chromosomeID, refAllele, altAllele;
  bool ignoreFamily;
  unsigned int ui;
//  char *stringArray = NULL;
  
  infile.seekg(8);
  infile.read((char *)&numSubjects, sizeof(int));
  infile.read((char *)&numSNPs, sizeof(int));
  infile.read((char *)&numGroups, sizeof(int));
  infile.read((char *)&subjectOptions, sizeof(unsigned int));
  ignoreFamily = (subjectOptions & 0x01) ? false : true;
  infile.read((char *)&snpOptions, sizeof(unsigned int));
  infile.read((char *)&uiStartSubjects, sizeof(unsigned int));
  infile.read((char *)&uiStartSNPs, sizeof(unsigned int));
  infile.read((char *)&uiStartDosage, sizeof(unsigned int));
  startSubjects = uiStartSubjects;
  startSNPs = uiStartSNPs;
  startDosage = uiStartDosage;
  infile.seekg(startSubjects);
  Rcpp::Rcout << "Subjects:\t" << numSubjects << "\tSNPs\t" << numSNPs << "\tGroups:\t" << numGroups << std::endl;
  Rcpp::Rcout << "Subject Options:\t" << std::hex << subjectOptions << "\tSNP Options:\t" << snpOptions << std::dec << std::endl;
  Rcpp::Rcout << "Start Subjects:\t" << uiStartSubjects << "\tStart SNPs:\t" << uiStartSNPs << "\tStart Dosage:\t" << uiStartDosage << std::endl;

  infile.read((char *)&subjectStringSize, sizeof(unsigned int));
  infile.read((char *)&familyStringSize, sizeof(unsigned int));
  Rcpp::Rcout << "Subject Array Size:\t" << subjectStringSize << "\tFamily Array Size:\t" << familyStringSize << std::endl;
  ReadStringArray(infile, subjectID, numSubjects, subjectStringSize);
  if (!ignoreFamily)
    ReadStringArray(infile, familyID, numSubjects, familyStringSize);
  subjects = Rcpp::List::create(Rcpp::Named("useFID") = ignoreFamily,
                                Rcpp::Named("FID") = familyID,
                                Rcpp::Named("IID") = subjectID);
  result["subjects"] = subjects;
  
  infile.read((char *)&snpStringSize, sizeof(unsigned int));
  infile.read((char *)&chromosomeStringSize, sizeof(unsigned int));
  infile.read((char *)&refAlleleStringSize, sizeof(unsigned int));
  infile.read((char *)&altAlleleStringSize, sizeof(unsigned int));
  Rcpp::Rcout << "SNP array size:\t" << snpStringSize << "\tChromosome array Size:\t" << chromosomeStringSize
              << "\tRef Allele array size:\t" << refAlleleStringSize << "\tAlt Allele array Size:\t" << altAlleleStringSize << std::endl;
  return true;
}
//' Function to produce summary of a binary dosage file used by GxEScan
//' 
//' Function to produce summary of a binary dosage file used by GxEScan
//' 
//' @param binaryDosageFilesname
//' Name of file with genetic data, normally ends with .bdosage
//' @param nSub
//' Number of subjects in the genetic data file. This is needed for formats 1, 2 and 3.
//' This value is returned for format 4 in the results
//' @param nSNPs
//' Number of SNPs in the genetic data file. This is needed for formats 1, 2 and 3.
//' This value is returned for format 4 in the results
//' @return
//' List with data needed by GxEScan
//' @export
// [[Rcpp::export]]
Rcpp::List GetBinaryDosageInformation(const std::string &binaryDosageFilename, const unsigned int nSub, const unsigned int nSNPs) {
  std::ifstream infile;
  const char magicWord[4] = { 'b', 'o', 's', 'e' };
  char readMagicWord[4];
  char readVersion[4];
  int format, version;
  unsigned int numSubjects, numSNPs;
  Rcpp::List result, subjects, snps;
  
  format = 0;
  version = 0;
  subjects = Rcpp::List::create(Rcpp::Named("useFID") = true,
                                Rcpp::Named("FID") = NULL,
                                Rcpp::Named("IID") = NULL);
  snps = Rcpp::List::create(Rcpp::Named("SNP") = NULL,
                            Rcpp::Named("Chromosome") = NULL,
                            Rcpp::Named("bp") = NULL,
                            Rcpp::Named("refAllele") = NULL,
                            Rcpp::Named("altAllele") = NULL,
                            Rcpp::Named("rsq") = NULL,
                            Rcpp::Named("callrate") = NULL);
  result = Rcpp::List::create(Rcpp::Named("success") = false,
                              Rcpp::Named("format") = format,
                              Rcpp::Named("version") = version,
                              Rcpp::Named("subjects") = subjects,
                              Rcpp::Named("SNPs") = snps,
                              Rcpp::Named("errorMessage") = "");
  
  infile.open(binaryDosageFilename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!infile.good()) {
    result["errorMessage"] = "Unable to open binary dosage file";
    return result;
  }
  infile.read(readMagicWord, 4);
  if (!infile.good()) {
    infile.close();
    result["errorMessage"] = "Error reading header of binary dosage file";
    return result;
  }
  if (std::memcmp(magicWord, readMagicWord, 4)) {
    infile.close();
    result["errorMessage"] = "File is not a binary dosage file";
    return result;
  }
  infile.read(readVersion, 4);
  if (!infile.good()) {
    infile.close();
    result["errorMessage"] = "Error reading format and version number";
    return result;
  }
  if (readVersion[0] || readVersion[2]) {
    infile.close();
    result["errorMessage"] = "Format and version number error";
    return result;
  }
  format = readVersion[1];
  version = readVersion[3];
  if (format < 1 || format > 4) {
    infile.close();
    result["errorMessage"] = "Unknown format";
    return result;
  }
  if (version < 1 || version > 2) {
    infile.close();
    result["errorMessage"] = "Unknown version";
    return result;
  }
  if (format > 2) {
    infile.read((char *)&numSubjects, sizeof(unsigned int));
    result["numSubjects"] = numSubjects;
    if (format > 3) {
      infile.read((char *)&numSNPs, sizeof(unsigned int));
      result["numSNPs"] = numSNPs;
      GetBinaryDosage4Info(infile, numSubjects, numSNPs, result);
    } else {
      CheckBinaryDosage3FileSize(infile, version, nSub, nSNPs);
    }
  } else {
    CheckBinaryDosageFileSize(infile, version, nSub, nSNPs);
  }
  
  result["success"] = true;
  result["format"] = format;
  result["version"] = version;
  infile.close();
  return result;
}
