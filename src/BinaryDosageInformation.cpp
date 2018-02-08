#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
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
//  Rcpp::Rcout << "Format:\tX." << version << '\t' << expectedSize << '\t' << actualSize << std::endl;
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
//    Rcpp::Rcout << "Format:\t3.1\t" << expectedSize << '\t' << actualSize << std::endl;
    if (expectedSize == actualSize)
      return true;
    Rcpp::Rcerr << "Binary dosage file is not of the expected size" << std::endl;
    return false;
  } else {
    infile.seekg(0, std::ios_base::end);
    actualSize = infile.tellg();
    expectedSize = 12 + sizeof(int) * nSNPs;
    infile.seekg(12);
    for (ui = 0; ui < nSNPs; ++ui) {
      infile.read((char *)&snpSize, sizeof(int));
//      Rcpp::Rcout << snpSize << std::endl;
      if (!infile.good()) {
        Rcpp::Rcerr << "Reached end of file:\t" << expectedSize << '\t' << actualSize << '\t' << snpSize << '\t' << infile.tellg() << std::endl;
        return false;
      }
      expectedSize += snpSize;
      infile.seekg(snpSize, std::ios_base::cur);
    }
    if (expectedSize != actualSize) {
      Rcpp::Rcerr << "Binary dosage file is not of the expected size\t" << expectedSize << '\t' << actualSize << std::endl;
      return false;
    }
//    Rcpp::Rcout << "Format:\t3.1\t" << expectedSize << '\t' << actualSize << std::endl;
  }
  
  return true;
}

void ReadStringArray(std::ifstream &infile, std::vector<std::string> &x, const unsigned int numStrings, const unsigned int arraySize) {
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

Rcpp::NumericMatrix ConvertToMatrix(const std::vector<double> &v, const unsigned int numSNPs, const unsigned int numGroups) {
  Rcpp::NumericMatrix res(numSNPs, numGroups);

  std::memcpy((char *)&res[0], (char *)&v[0], numSNPs * numGroups * sizeof(double));
  return res;
}

bool GetBinaryDosage4Info(std::ifstream &infile, unsigned int &numSubjects, unsigned int &numSNPs, Rcpp::List &result) {
  Rcpp::List subjects;
  Rcpp::List snps;
  Rcpp::DataFrame info;
  unsigned int numGroups;
  unsigned int subjectOptions;
  unsigned int snpOptions;
  unsigned int uiStartSubjects, uiStartSNPs, uiStartDosage;
  std::streampos startSubjects, startSNPs, startDosage;
  unsigned int subjectStringSize, familyStringSize, snpStringSize, chromosomeStringSize, refAlleleStringSize, altAlleleStringSize;
  std::vector<std::string> subjectID, familyID, snpID, chromosomeID, refAllele, altAllele;
  std::vector<int> bp;
  std::vector<double> aaf, maf, avgCall, rsq, tempVec;
  std::vector<int> groupSize;
  Rcpp::NumericMatrix aaf2, maf2, avgCall2, rsq2;
  bool ignoreFamily;
  unsigned int ui, uj, uk;
  unsigned int readStart;

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
  groupSize.resize(numGroups);
  infile.read((char *)&groupSize[0], numGroups*sizeof(int));
  infile.seekg(startSubjects);
//  Rcpp::Rcout << "Subjects:\t" << numSubjects << "\tSNPs\t" << numSNPs << "\tGroups:\t" << numGroups << std::endl;
//  Rcpp::Rcout << "Group Sizes:";
//  for (ui = 0; ui < numGroups; ++ui)
//    Rcpp::Rcout << '\t' << groupSize[ui];
//  Rcpp::Rcout << std::endl;
//  Rcpp::Rcout << "Subject Options:\t" << std::hex << subjectOptions << "\tSNP Options:\t" << snpOptions << std::dec << std::endl;
//  Rcpp::Rcout << "Start Subjects:\t" << uiStartSubjects << "\tStart SNPs:\t" << uiStartSNPs << "\tStart Dosage:\t" << uiStartDosage << std::endl;

  infile.read((char *)&subjectStringSize, sizeof(unsigned int));
  infile.read((char *)&familyStringSize, sizeof(unsigned int));
//  Rcpp::Rcout << "Subject Array Size:\t" << subjectStringSize << "\tFamily Array Size:\t" << familyStringSize << std::endl;
  ReadStringArray(infile, subjectID, numSubjects, subjectStringSize);
  if (!ignoreFamily)
    ReadStringArray(infile, familyID, numSubjects, familyStringSize);
  subjects = Rcpp::List::create(Rcpp::Named("useFID") = !ignoreFamily,
                                Rcpp::Named("FID") = familyID,
                                Rcpp::Named("IID") = subjectID);
//                                Rcpp::Named("stringsAsFactors") = false);
  result["subjects"] = subjects;
  
  infile.read((char *)&snpStringSize, sizeof(unsigned int));
  infile.read((char *)&chromosomeStringSize, sizeof(unsigned int));
  infile.read((char *)&refAlleleStringSize, sizeof(unsigned int));
  infile.read((char *)&altAlleleStringSize, sizeof(unsigned int));
//  Rcpp::Rcout << "SNP array size:\t" << snpStringSize << "\tChromosome array Size:\t" << chromosomeStringSize
//              << "\tRef Allele array size:\t" << refAlleleStringSize << "\tAlt Allele array Size:\t" << altAlleleStringSize << std::endl;
  if (snpOptions & 0x02) 
    ReadStringArray(infile, snpID, numSNPs, snpStringSize);
  if (snpOptions & 0x08) {
    ReadStringArray(infile, chromosomeID, 1, chromosomeStringSize);
    chromosomeID.resize(numSNPs);
    for (ui = 1; ui < numSNPs; ++ui)
      chromosomeID[ui] = chromosomeID[0];
  } else {
    ReadStringArray(infile, chromosomeID, numSNPs, chromosomeStringSize);
  }
  bp.resize(numSNPs);
  infile.read((char *)&bp[0], numSNPs * sizeof(int));
  if ((snpOptions & 0x02) == 0) {
    snpID.resize(numSNPs);
    for (ui = 0; ui < numSNPs; ++ui)
      snpID[ui] = chromosomeID[ui] + ":" + std::to_string(bp[ui]);
  }
  if (snpOptions &0x0020)
    ReadStringArray(infile, refAllele, numSNPs, refAlleleStringSize);
  if (snpOptions &0x0040)
    ReadStringArray(infile, altAllele, numSNPs, altAlleleStringSize);
  
  if (snpOptions &0x0080)
    aaf.resize(numSNPs * numGroups);
  if (snpOptions &0x0100)
    maf.resize(numSNPs * numGroups);
  if (snpOptions &0x0200)
    avgCall.resize(numSNPs * numGroups);
  if (snpOptions &0x0400)
    rsq.resize(numSNPs * numGroups);
  tempVec.resize(numSNPs * numGroups);
  
  if (snpOptions &0x0080) {
    infile.read((char *)&tempVec[0], numGroups * numSNPs * sizeof(double));
    
    readStart = 0;
    for (ui = 0; ui < numGroups; ++ui, readStart += numSNPs) {
      uj = ui;
      for (uk = 0; uk < numSNPs; ++uk, uj += numGroups)
        aaf[readStart + uk] = tempVec[uj];
    }
    aaf2 = ConvertToMatrix(aaf, numSNPs, numGroups);
  }
  if (snpOptions &0x0100) {
    infile.read((char *)&tempVec[0], numGroups * numSNPs * sizeof(double));
    readStart = 0;
    for (ui = 0; ui < numGroups; ++ui, readStart += numSNPs) {
      uj = ui;
      for (uk = 0; uk < numSNPs; ++uk, uj += numGroups)
        maf[readStart + uk] = tempVec[uj];
    }
    maf2 = ConvertToMatrix(maf, numSNPs, numGroups);
  }
  if (snpOptions &0x0200) {
    infile.read((char *)&tempVec[0], numGroups * numSNPs * sizeof(double));
    readStart = 0;
    for (ui = 0; ui < numGroups; ++ui, readStart += numSNPs) {
      uj = ui;
      for (uk = 0; uk < numSNPs; ++uk, uj += numGroups)
        avgCall[readStart + uk] = tempVec[uj];
    }
    avgCall2 = ConvertToMatrix(avgCall, numSNPs, numGroups);
  }
  if (snpOptions &0x0400) {
    infile.read((char *)&tempVec[0], numGroups * numSNPs * sizeof(double));
    readStart = 0;
    for (ui = 0; ui < numGroups; ++ui, readStart += numSNPs) {
      uj = ui;
      for (uk = 0; uk < numSNPs; ++uk, uj += numGroups)
        rsq[readStart + uk] = tempVec[uj];
    }
    rsq2 = ConvertToMatrix(rsq, numSNPs, numGroups);
  }
  
  info = DataFrame::create(Rcpp::Named("ID") = snpID,
                           Rcpp::Named("Chromosome") = chromosomeID,
                           Rcpp::Named("bp") = bp,
                           Rcpp::Named("refAllele") = refAllele,
                           Rcpp::Named("altAllele") = altAllele,
                           Rcpp::Named("stringsAsFactors") = false);
//  Rcpp::Rcout << "MAF size:\t" << maf.size() << std::endl;
//  Rcpp::Rcout << infile.tellg() << std::endl;
  snps = Rcpp::List::create(Rcpp::Named("Info") = info,
                            Rcpp::Named("aaf") = aaf2,
                            Rcpp::Named("maf") = maf2,
                            Rcpp::Named("avgCall") = avgCall2,
                            Rcpp::Named("rsq") = rsq2);
  result["SNPs"] = snps;
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
  Rcpp::List result;//, subjects, snps;
  
  format = 0;
  version = 0;
/*
  subjects = Rcpp::List::create(Rcpp::Named("useFID") = true,
                                Rcpp::Named("FID") = NULL,
                                Rcpp::Named("IID") = NULL);
  snps = Rcpp::List::create(Rcpp::Named("SNP") = NULL,
                            Rcpp::Named("Chromosome") = NULL,
                            Rcpp::Named("bp") = NULL,
                            Rcpp::Named("refAllele") = NULL,
                            Rcpp::Named("altAllele") = NULL,
                            Rcpp::Named("rsq") = NULL,
                            Rcpp::Named("maf") = NULL);
 */
  result = Rcpp::List::create(Rcpp::Named("geneticFile") = binaryDosageFilename,
                              Rcpp::Named("success") = false,
                              Rcpp::Named("format") = format,
                              Rcpp::Named("version") = version,
//                              Rcpp::Named("subjects") = subjects,
//                              Rcpp::Named("SNPs") = snps,
                              Rcpp::Named("errorMessage") = "",
                              Rcpp::Named("numSubjects") = nSub,
                              Rcpp::Named("numSNPs") = nSNPs);
  
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
