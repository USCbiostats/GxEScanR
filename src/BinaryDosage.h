#ifndef BINARYDOSAGE_H
#define BINARYDOSAGE_H 1

#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

class CBinaryDosage : public CGeneticData {
protected:
  std::ifstream m_infile;
  std::string m_geneticFilename;
  int m_version;
  int m_subversion;
  
  CBinaryDosage();
  CBinaryDosage(const int _numSubjects, const int _numSNPs, bool _geneticProbabilities, std::string &_geneticFilename, const int _version, const int _subversion);
public:
  virtual ~CBinaryDosage() { m_infile.close(); }
};

class CBinaryDosageFormat1_1 : public CBinaryDosage {
public:
  CBinaryDosageFormat1_1(std::string &_geneticFilename, const int numSubjects, const int numSNPs);
  virtual ~CBinaryDosageFormat1_1() {}

  virtual int GetFirst();
  virtual int GetNext() { return 1; }
  virtual int GetSNP(unsigned int n) { return 1; }
};

/*
class CBinaryDosageFormat4_2 : public CBinaryDosage {
protected:
};
*/

#endif