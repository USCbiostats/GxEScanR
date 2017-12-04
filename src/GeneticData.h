#ifndef GENETICDATA_H
#define GENETIDDATA_H 1

int TestGeneticData(Rcpp::List &_geneticData);

class CGeneticData {
protected:
  std::string m_errorMessage;
  
  CGeneticData() { m_errorMessage = ""; }
  virtual ~CGeneticData() {}
public:
  virtual int CheckValidity() = 0;
  const std::string &ErrorMessage() { return m_errorMessage; }
};
#endif