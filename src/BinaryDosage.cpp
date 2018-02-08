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
}

int CBinaryDosageFormat1_1::GetFirst() {
  if (!m_valid)
    return 1;
  m_infile.clear();
  m_infile.seekg(8);
//  m_infile.read();
  return 0;
}