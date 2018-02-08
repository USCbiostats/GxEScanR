#include <string>
#include <RcppArmadillo.h>
#include "GeneticData.h"
#include "PlinkBinary.h"
#include "BinaryDosage.h"

CGeneticData::CGeneticData(const int _numSubjects, const int _numSNPs, bool _measured, bool _geneticProbabilities) {
  m_errorMessage = "";
  m_numSubjects = _numSubjects;
  m_numSNPs = _numSNPs;
  m_valid = false;
  m_bMeasured = _measured;
  m_bGeneticProbabilities = _geneticProbabilities;
  m_dosages.resize(m_numSubjects);
  if (m_bGeneticProbabilities)
    m_probabilities.resize(m_numSubjects, 3);
}

int CGeneticData::GetSNP(unsigned int n) {
  if (!m_valid) {
    m_errorMessage = "Data for genetic file is invalid";
    return 1;
  }
  if (n == 0 || n > m_numSNPs) {
    m_errorMessage = "Invalid SNP value";
    return 1;
  }
  return 0;
}