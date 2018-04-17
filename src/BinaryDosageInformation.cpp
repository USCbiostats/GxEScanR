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
  Rcpp::stop("Binary dosage file is not of the expected size");
//  Rcpp::Rcerr << "Binary dosage file is not of the expected size" << std::endl;
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
    Rcpp::stop("Binary dosage file is not of the expected size");
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
        Rcpp::stop("Binary dosage file is not of the expected size");
        //        Rcpp::Rcerr << "Reached end of file:\t" << expectedSize << '\t' << actualSize << '\t' << snpSize << '\t' << infile.tellg() << std::endl;
        return false;
      }
      expectedSize += snpSize;
      infile.seekg(snpSize, std::ios_base::cur);
    }
    if (expectedSize != actualSize) {
      Rcpp::stop("Binary dosage file is not of the expected size");
      //      Rcpp::Rcerr << "Binary dosage file is not of the expected size\t" << expectedSize << '\t' << actualSize << std::endl;
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
  Rcpp::DataFrame subjects;
  Rcpp::DataFrame snps;
  Rcpp::DataFrame SNPinfo;
  unsigned int numGroups;
  unsigned int subjectOptions;
  unsigned int snpOptions;
  unsigned int uiStartSubjects, uiStartSNPs, uiStartDosage;
  std::streampos startSubjects, startSNPs, startDosage;
  unsigned int subjectStringSize, familyStringSize, snpStringSize, chromosomeStringSize, refAlleleStringSize, altAlleleStringSize;
  std::vector<std::string> familyID, subjectID, snpID, chromosomeID, refAllele, altAllele;
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

  // Read in the family and subject IDs
  infile.read((char *)&subjectStringSize, sizeof(unsigned int));
  infile.read((char *)&familyStringSize, sizeof(unsigned int));
//  Rcpp::Rcout << "Subject Array Size:\t" << subjectStringSize << "\tFamily Array Size:\t" << familyStringSize << std::endl;
  ReadStringArray(infile, subjectID, numSubjects, subjectStringSize);
  if (!ignoreFamily)
    ReadStringArray(infile, familyID, numSubjects, familyStringSize);
  else
    familyID.resize(numSubjects);
  subjects = Rcpp::DataFrame::create(Rcpp::Named("FID") = familyID,
                                    Rcpp::Named("IID") = subjectID,
                                    Rcpp::Named("stringsAsFactors") = false);
  result["numSubjects"] = numSubjects;
  result["subjects"] = subjects;

  // Read in data about the SNPs  
  infile.read((char *)&snpStringSize, sizeof(unsigned int));
  infile.read((char *)&chromosomeStringSize, sizeof(unsigned int));
  infile.read((char *)&refAlleleStringSize, sizeof(unsigned int));
  infile.read((char *)&altAlleleStringSize, sizeof(unsigned int));
//  Rcpp::Rcout << "SNP array size:\t" << snpStringSize << "\tChromosome array Size:\t" << chromosomeStringSize
//              << "\tRef Allele array size:\t" << refAlleleStringSize << "\tAlt Allele array Size:\t" << altAlleleStringSize << std::endl;

  // Read in the SNP names - SNP names will be set to chromosome:location if not read in here
  if (snpOptions & 0x02) 
    ReadStringArray(infile, snpID, numSNPs, snpStringSize);
  else
    snpID.resize(numSNPs);

  // Read in the chromosome - Why is this required? It isn't if SNP names are provided. 
  if (snpOptions & 0x08) {
    ReadStringArray(infile, chromosomeID, 1, chromosomeStringSize);
    // Resize shouldn't change initial value, I hope ???
    chromosomeID.resize(numSNPs);
    for (ui = 1; ui < numSNPs; ++ui)
      chromosomeID[ui] = chromosomeID[0];
  } else {
    ReadStringArray(infile, chromosomeID, numSNPs, chromosomeStringSize);
  }
  
  // Read in the SNP locations - again it isn't required when SNP names are provided.
  bp.resize(numSNPs);
  infile.read((char *)&bp[0], numSNPs * sizeof(int));
  
  // Assign SNP names if they weren't read earlier
  if ((snpOptions & 0x02) == 0) {
    snpID.resize(numSNPs);
    for (ui = 0; ui < numSNPs; ++ui)
      snpID[ui] = chromosomeID[ui] + ":" + std::to_string(bp[ui]);
  }
  
  // Read in reference and alternate alleles if provided
  if (snpOptions &0x0020)
    ReadStringArray(infile, refAllele, numSNPs, refAlleleStringSize);
  else
    refAllele.resize(numSNPs);
  if (snpOptions &0x0040)
    ReadStringArray(infile, altAllele, numSNPs, altAlleleStringSize);
  else
    altAllele.resize(numSNPs);
  
  // Read in imputation info about SNPs
  aaf.resize(numSNPs * numGroups);
  maf.resize(numSNPs * numGroups);
  avgCall.resize(numSNPs * numGroups);
  rsq.resize(numSNPs * numGroups);
  tempVec.resize(numSNPs * numGroups);

  // The follow sections of code are required because the information was
  // saved in a row dominant format
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
  
  snps = DataFrame::create(Rcpp::Named("SNP") = snpID,
                           Rcpp::Named("CHR") = chromosomeID,
                           Rcpp::Named("BP") = bp,
                           Rcpp::Named("A1") = refAllele,
                           Rcpp::Named("A2") = altAllele,
                           Rcpp::Named("stringsAsFactors") = false);
//  Rcpp::Rcout << "MAF size:\t" << maf.size() << std::endl;
//  Rcpp::Rcout << infile.tellg() << std::endl;
  SNPinfo = Rcpp::DataFrame::create(Rcpp::Named("aaf") = aaf2,
                                    Rcpp::Named("maf") = maf2,
                                    Rcpp::Named("avgCall") = avgCall2,
                                    Rcpp::Named("rsq") = rsq2);
  result["snps"] = snps;
  result["SNPinfo"] = SNPinfo;
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

  Rcpp::DataFrame subjects = Rcpp::DataFrame::create(Rcpp::Named("FID") = Rcpp::CharacterVector(),
                                     Rcpp::Named("IID") = Rcpp::CharacterVector(),
                                     Rcpp::Named("stringsAsFactors") = false);
  Rcpp::DataFrame snps = Rcpp::DataFrame::create(Rcpp::Named("SNP") = Rcpp::CharacterVector(),
                                                 Rcpp::Named("CHR") = Rcpp::CharacterVector(),
                                                 Rcpp::Named("BP") = Rcpp::IntegerVector(),
                                                 Rcpp::Named("A1") = Rcpp::CharacterVector(),
                                                 Rcpp::Named("A2") = Rcpp::CharacterVector(),
                                                 Rcpp::Named("rsq") = Rcpp::NumericVector(),
                                                 Rcpp::Named("maf") = Rcpp::NumericVector(),
                                                 Rcpp::Named("stringsAsFactors") = false);
  Rcpp::DataFrame SNPInfo = Rcpp::DataFrame::create(Rcpp::Named("aaf") = Rcpp::NumericVector(),
                                                    Rcpp::Named("maf") = Rcpp::NumericVector(),
                                                    Rcpp::Named("avgCall") = Rcpp::NumericVector(),
                                                    Rcpp::Named("rsq") = Rcpp::NumericVector());

  result = Rcpp::List::create(Rcpp::Named("filetype") = "BinaryDosage",
                              Rcpp::Named("filename") = binaryDosageFilename,
                              Rcpp::Named("format") = format,
                              Rcpp::Named("version") = version,
                              Rcpp::Named("numSubjects") = 0,
                              Rcpp::Named("usesFID") = true,
                              Rcpp::Named("subjects") = subjects,
                              Rcpp::Named("numSNPs") = 0,
                              Rcpp::Named("snps") = snps,
                              Rcpp::Named("SNPinfo") = SNPInfo);

  // Open the file  
  infile.open(binaryDosageFilename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!infile.good())
    Rcpp::stop("Unable to open binary dosage file");
  // Check the header
  // Check magic word
  infile.read(readMagicWord, 4);
  if (!infile.good()) {
    infile.close();
    Rcpp::stop("Error reading header of binary dosage file");
  }
  if (std::memcmp(magicWord, readMagicWord, 4)) {
    infile.close();
    Rcpp::stop("File is not a binary dosage file");
  }
  
  // Check version
  infile.read(readVersion, 4);
  if (!infile.good()) {
    infile.close();
    Rcpp::stop("Error reading format and version number");
  }
  if (readVersion[0] || readVersion[2]) {
    infile.close();
    Rcpp::stop("Format and version number error");
  }
  format = readVersion[1];
  version = readVersion[3];
  if (format < 1 || format > 4) {
    infile.close();
    Rcpp::stop("Unknown format");
  }
  if (version < 1 || version > 2) {
    infile.close();
    Rcpp::stop("Unknown version");
  }
  result["format"] = format;
  result["version"] = version;
  
  // If format is greater than 2, read the number of subjects
  if (format > 2) {
    infile.read((char *)&numSubjects, sizeof(unsigned int));
    result["numSubjects"] = numSubjects;
    // If format is greater than 3, read the number of SNPs
    // and the subject and SNP data
    if (format > 3) {
      infile.read((char *)&numSNPs, sizeof(unsigned int));
      result["numSNPs"] = numSNPs;
      GetBinaryDosage4Info(infile, numSubjects, numSNPs, result);
    } else {
      // Check the file size for format 3 - different than 1 and 2
      CheckBinaryDosage3FileSize(infile, version, nSub, nSNPs);
    }
  } else {
    // Check the file size for formats 1 and 2
    CheckBinaryDosageFileSize(infile, version, nSub, nSNPs);
  }
  
  infile.close();
  return result;
}
