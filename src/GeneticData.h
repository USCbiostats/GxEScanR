#ifndef GENETICDATA_H
#define GENETIDDATA_H 1

int TestGeneticData(Rcpp::List &_geneticData);

class CGeneticData {
protected:
  std::string m_errorMessage;
  unsigned int m_numSubjects;
  unsigned int m_numSNPs;
  
  CGeneticData();
  virtual ~CGeneticData() {}
public:
  virtual int CheckValidity() = 0;
  const std::string &ErrorMessage() { return m_errorMessage; }
  
  unsigned int NumSubjects() const { return m_numSubjects; }
  unsigned int NumSNPs() const { return m_numSNPs; }
};
#endif