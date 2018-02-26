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
  CBinaryDosage(const int _numSubjects, const int _numSNPs, bool _geneticProbabilities, const std::string &_geneticFilename, const int _version, const int _subversion);
public:
  virtual ~CBinaryDosage() { m_infile.close(); }
};

class CBinaryDosageFormat1_1 : public CBinaryDosage {
protected:
  std::vector<unsigned short> m_readBuffer;
  virtual void AssignDosages();
public:
  CBinaryDosageFormat1_1(const std::string &_geneticFilename, const int numSubjects, const int numSNPs);
  virtual ~CBinaryDosageFormat1_1() {}

  virtual int GetFirst();
  virtual int GetNext();
  virtual int GetSNP(unsigned int n);
};

class CBinaryDosageFormat1_2 : public CBinaryDosage {
protected:
  std::vector<unsigned short> m_readBuffer;
  virtual void AssignDosages();
public:
  CBinaryDosageFormat1_2(const std::string &_geneticFilename, const int numSubjects, const int numSNPs);
  virtual ~CBinaryDosageFormat1_2() {}
  
  virtual int GetFirst();
  virtual int GetNext();
  virtual int GetSNP(unsigned int n);
};

class CBinaryDosageFormat2_1 : public CBinaryDosage {
protected:
  std::vector<unsigned short> m_readBuffer;
  virtual void AssignDosages();
public:
  CBinaryDosageFormat2_1(const std::string &_geneticFilename, const int numSubjects, const int numSNPs);
  virtual ~CBinaryDosageFormat2_1() {}
  
  virtual int GetFirst();
  virtual int GetNext();
  virtual int GetSNP(unsigned int n);
};

class CBinaryDosageFormat2_2 : public CBinaryDosage {
protected:
  std::vector<unsigned short> m_readBuffer;
  virtual void AssignDosages();
public:
  CBinaryDosageFormat2_2(const std::string &_geneticFilename, const int numSubjects, const int numSNPs);
  virtual ~CBinaryDosageFormat2_2() {}
  
  virtual int GetFirst();
  virtual int GetNext();
  virtual int GetSNP(unsigned int n);
};

class CBinaryDosageFormat3_1 : public CBinaryDosage {
protected:
  std::vector<unsigned short> m_readBuffer;
  virtual void AssignDosages();
public:
  CBinaryDosageFormat3_1(const std::string &_geneticFilename, const int numSubjects, const int numSNPs);
  virtual ~CBinaryDosageFormat3_1() {}
  
  virtual int GetFirst();
  virtual int GetNext();
  virtual int GetSNP(unsigned int n);
};

class CBinaryDosageFormat3_2 : public CBinaryDosage {
protected:
  std::vector<unsigned short> m_readBuffer;
  virtual void AssignDosages();
public:
  CBinaryDosageFormat3_2(const std::string &_geneticFilename, const int numSubjects, const int numSNPs);
  virtual ~CBinaryDosageFormat3_2() {}
  
  virtual int GetFirst();
  virtual int GetNext();
  virtual int GetSNP(unsigned int n);
};

/*
class CBinaryDosageFormat4_2 : public CBinaryDosage {
protected:
};
*/

#endif