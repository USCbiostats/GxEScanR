#ifndef BINARYDOSAGE_H
#define BINARYDOSAGE_H 1

#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

Rcpp::List BinaryDosageInfo(std::string &geneticFile, std::string &mapFile, std::string &familyFile);

int GetBinaryDosageFormat(std::string &binaryDosageFile, int &format, int &version, std::string &errorMessage);

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
};

class CBinaryDosageFormat1 : public CBinaryDosage {
protected:
  std::string m_mapFile;
  std::string m_familyFile;
public:
  CBinaryDosageFormat1(std::string &_geneticFile, std::string &_mapFile, std::string &_familyFile);
  CBinaryDosageFormat1(Rcpp::List &_binaryDosageInfo);
  virtual ~CBinaryDosageFormat1() {}
};

class CBinaryDosageFormat1_1 : public CBinaryDosageFormat1 {
public:
  CBinaryDosageFormat1_1(std::string &_geneticFile, std::string &_mapFile, std::string &_familyFile);
  CBinaryDosageFormat1_1(Rcpp::List &_binaryDosageInfo);
  virtual ~CBinaryDosageFormat1_1() {}
};

#endif