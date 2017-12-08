#ifndef GENETICDATA_H
#define GENETICDATA_H 1

class CGeneticData;

CGeneticData *TestGeneticData(Rcpp::List &_geneticData);

class CGeneticData {
protected:
  std::string m_errorMessage;
  unsigned int m_numSubjects;
  unsigned int m_numSNPs;
  
  CGeneticData();
public:
  virtual ~CGeneticData() {}

  virtual int CheckValidity() = 0;

  const std::string &ErrorMessage() { return m_errorMessage; }
  
  unsigned int NumSubjects() const { return m_numSubjects; }
  unsigned int NumSNPs() const { return m_numSNPs; }
};

#endif