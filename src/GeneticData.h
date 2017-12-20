#ifndef GENETICDATA_H
#define GENETICDATA_H 1

class CGeneticData;

CGeneticData *TestGeneticData(Rcpp::List &_geneticData);

class CGeneticData {
protected:
  std::string m_errorMessage;
  unsigned int m_numSubjects;
  unsigned int m_numSNPs;
  bool m_valid;
  
  CGeneticData();
public:
  virtual ~CGeneticData() {}

  virtual int CheckValidity() = 0;

  const std::string &ErrorMessage() { return m_errorMessage; }
  
  unsigned int NumSubjects() const { return m_numSubjects; }
  unsigned int NumSNPs() const { return m_numSNPs; }
  bool Valid() const { return m_valid; }
  
  virtual int GetFirst() = 0;
  virtual int GetNext() = 0;
  virtual int GetSNP(unsigned int n);
};

#endif