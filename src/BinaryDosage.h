#ifndef BINARYDOSAGE_H
#define BINARYDOSAGE_H 1

#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

Rcpp::List BinaryDosageInfo(const std::string &geneticFilename, const std::string &mapFilename, const std::string &familyFilename);

int GetBinaryDosageFormat(const std::string &binaryDosageFilename, int &format, int &version, std::string &errorMessage);

int GetNumberOfSubjectsAndSNPs(const std::string &familyFilename, const std::string &mapFilename, unsigned int &numSubjects, unsigned int &numSNPs, std::string &errorMessage);

int GetNumberOfSubjectsAndSNPsV4(const std::string &geneticFilename, unsigned int &numSubjects, unsigned int &numSNPs, unsigned int &groups, std::string &errorMessage);

int CheckBinaryDosageFileSize(const std::string &geneticFilename, int version, unsigned int numSubjects, unsigned int numSNPs, std::string &errorMessage);

class CBinaryDosage : public CGeneticData {
protected:
  std::ifstream m_infile;
  std::string m_geneticFilename;
  int m_version;
  int m_subversion;
  
  CBinaryDosage();
  CBinaryDosage(std::string &_geneticFilename);
public:
  virtual ~CBinaryDosage() { m_infile.close(); }

  virtual int GetFirst() { return 1; }
  virtual int GetNext() { return 1; }
};

class CBinaryDosageFormat1 : public CBinaryDosage {
protected:
  std::string m_mapFilename;
  std::string m_familyFilename;
public:
  CBinaryDosageFormat1(std::string &_geneticFilename, std::string &_mapFilename, std::string &_familyFilename);
  CBinaryDosageFormat1(Rcpp::List &_binaryDosageInfo);
  virtual ~CBinaryDosageFormat1() {}

  virtual int CheckVersion() = 0;
  virtual int CheckValidity();
};

class CBinaryDosageFormat1_1 : public CBinaryDosageFormat1 {
public:
  CBinaryDosageFormat1_1(std::string &_geneticFilename, std::string &_mapFilename, std::string &_familyFilename);
  CBinaryDosageFormat1_1(Rcpp::List &_binaryDosageInfo);
  virtual ~CBinaryDosageFormat1_1() {}
  
  virtual int CheckVersion();
};
/*
class CBinaryDosageFormat4_2 : public CBinaryDosage {
protected:
};
*/

#endif