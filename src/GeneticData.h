#ifndef GENETICDATA_H
#define GENETICDATA_H 1

#include <Rcpp.h>

class CGeneticData {
protected:
  std::string m_errorMessage;
  unsigned int m_numSubjects;
  unsigned int m_numSNPs;
  bool m_valid;
  bool m_bMeasured;
  bool m_bGeneticProbabilities;
  double *m_doseData;
  double *m_dosages;
  double *m_probabilities[3];
  
  CGeneticData(const int _numSubjects, const int _numSNPs, bool _measured, bool _geneticProbabilities);
  CGeneticData() {}
public:
  virtual ~CGeneticData();

  const std::string &ErrorMessage() { return m_errorMessage; }
  unsigned int NumSubjects() const { return m_numSubjects; }
  unsigned int NumSNPs() const { return m_numSNPs; }
  bool Valid() const { return m_valid; }
  bool Measured() const { return m_bMeasured; }
  bool GeneticProbabilities() const { return m_bGeneticProbabilities; }
  const double *Dosages() const { return m_dosages; }
  const double *Probabilities(int n) const { return m_probabilities[n]; }
  
  virtual int GetFirst() = 0;
  virtual int GetNext() = 0;
  virtual int GetSNP(unsigned int n);
};

#endif