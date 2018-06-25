#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
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
  if (m_bGeneticProbabilities) {
    m_doseData = new double[4 * m_numSubjects];
    m_dosages = m_doseData;
    m_probabilities[0] = m_dosages + m_numSubjects;
    m_probabilities[1] = m_probabilities[0] + m_numSubjects;
    m_probabilities[2] = m_probabilities[1] + m_numSubjects;
  } else {
    m_doseData = new double[m_numSubjects];
    m_dosages = m_doseData;
    m_probabilities[0] = NULL;
    m_probabilities[1] = NULL;
    m_probabilities[2] = NULL;
  }
}

CGeneticData::~CGeneticData() {
  if (m_doseData)
    delete [] m_doseData;
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