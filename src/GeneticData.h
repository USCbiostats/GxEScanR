#ifndef GENETICDATA_H
#define GENETICDATA_H 1

#include <RcppArmadillo.h>

class CGeneticData {
protected:
  std::string m_errorMessage;
  unsigned int m_numSubjects;
  unsigned int m_numSNPs;
  bool m_valid;
  bool m_bMeasured;
  bool m_bGeneticProbabilities;
  arma::Col<double> m_dosages;
  arma::Mat<double> m_probabilities;
  
  CGeneticData(const int _numSubjects, const int _numSNPs, bool _measured, bool _geneticProbabilities);
  CGeneticData() {}
public:
  virtual ~CGeneticData() {}

  const std::string &ErrorMessage() { return m_errorMessage; }
  unsigned int NumSubjects() const { return m_numSubjects; }
  unsigned int NumSNPs() const { return m_numSNPs; }
  bool Valid() const { return m_valid; }
  bool Measured() const { return m_bMeasured; }
  bool GeneticProbabilities() const { return m_bGeneticProbabilities; }
  const arma::Col<double> Dosages() const { return m_dosages; }
  const arma::Mat<double> Probabilities() const { return m_probabilities; }
  
  virtual int GetFirst() = 0;
  virtual int GetNext() = 0;
  virtual int GetSNP(unsigned int n);
};

#endif