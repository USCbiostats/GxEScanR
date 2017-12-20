#ifndef BINARYDOSAGE_H
#define BINARYDOSAGE_H 1

#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

Rcpp::List BinaryDosageInfo(const std::string &geneticFile, const std::string &mapFile, const std::string &familyFile);

int GetBinaryDosageFormat(const std::string &binaryDosageFile, int &format, int &version, std::string &errorMessage);

int GetNumberOfSubjectsAndSNPs(const std::string &familyFile, const std::string &mapFile, int &numSubjects, int &numSNPs, std::string &errorMessage);

int GetNumberOfSubjectsAndSNPsV4(const std::string &geneticFile, int &numSubjects, int &numSNPs, int &groups, std::string &errorMessage);

int CheckBinaryDosageFileSize(const std::string &geneticFile, int version, int numSubjects, int numSNPs, std::string &errorMessage);

class CBinaryDosage : public CGeneticData {
protected:
  std::ifstream m_infile;
  std::string m_geneticFile;
  int m_version;
  int m_subversion;
  
  CBinaryDosage();
public:
  CBinaryDosage(std::string &_geneticFile);
  virtual ~CBinaryDosage() { m_infile.close(); }

  virtual int GetFirst() { return 1; }
  virtual int GetNext() { return 1; }
};

class CBinaryDosageFormat1 : public CBinaryDosage {
protected:
  std::string m_mapFile;
  std::string m_familyFile;
public:
  CBinaryDosageFormat1(std::string &_geneticFile, std::string &_mapFile, std::string &_familyFile);
  CBinaryDosageFormat1(Rcpp::List &_binaryDosageInfo);
  virtual ~CBinaryDosageFormat1() {}

  virtual int CheckVersion() = 0;
  virtual int CheckValidity();
};

class CBinaryDosageFormat1_1 : public CBinaryDosageFormat1 {
public:
  CBinaryDosageFormat1_1(std::string &_geneticFile, std::string &_mapFile, std::string &_familyFile);
  CBinaryDosageFormat1_1(Rcpp::List &_binaryDosageInfo);
  virtual ~CBinaryDosageFormat1_1() {}
  
  virtual int CheckVersion();
};

#endif