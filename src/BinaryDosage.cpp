#include <istream>
#include <fstream>
#include <string>
#include <RcppArmadillo.h>
#include "BinaryDosage.h"

CBinaryDosage::CBinaryDosage(const int _numSubjects, const int _numSNPs, bool _geneticProbabilities, std::string &_geneticFilename, const int _version, const int _subversion)
  : CGeneticData(_numSubjects, _numSNPs, false, _geneticProbabilities) {
  m_geneticFilename = _geneticFilename;
  m_infile.open(_geneticFilename.c_str(), std::ios_base::in | std::ios_base::binary);
}

CBinaryDosageFormat1_1::CBinaryDosageFormat1_1(std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, false, _geneticFilename, 1, 1) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x1, 0x0, 0x1};
  char readHeader[8];
  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
    if (std::memcmp(header, readHeader, 8))
      m_valid = true;
    else
      m_errorMessage = "File is not a binary dosage format 1.1 file";
  }
  m_readBuffer.resize(m_numSubjects);
}

void CBinaryDosageFormat1_1::AssignDosages() {
  unsigned short *u;
  double *d;
  unsigned int ui;
  
  for (ui = 0, u = &m_readBuffer[0], d = &m_dosages[0]; ui < m_numSubjects; ++ui, ++u, ++d) {
    if (*u == 0xffff)
      *d = NA_REAL;
    else
      *d = *u / (double)0xfffe;
  }
}

int CBinaryDosageFormat1_1::GetFirst() {
  return GetSNP(1);
}

int CBinaryDosageFormat1_1::GetNext() {
  if (!m_valid)
    return 1;
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}

int CBinaryDosageFormat1_1::GetSNP(unsigned int n) {
  if (CGeneticData::GetSNP(n))
    return 1;
  m_infile.clear();
  m_infile.seekg(8 + (n - 1) * m_numSNPs * sizeof(unsigned int));
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}