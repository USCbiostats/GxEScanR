#include <istream>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include "BinaryDosage.h"

CBinaryDosage::CBinaryDosage(const int _numSubjects, const int _numSNPs, bool _geneticProbabilities, const std::string &_geneticFilename, const int _version, const int _subversion)
  : CGeneticData(_numSubjects, _numSNPs, false, _geneticProbabilities) {
  m_geneticFilename = _geneticFilename;
  m_infile.open(_geneticFilename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!m_infile.good())
    Rcpp::Rcerr << "Unable to open file:\t" << _geneticFilename << std::endl;
}

// ******************   CBinaryDosageFormat1_1   **********************/

CBinaryDosageFormat1_1::CBinaryDosageFormat1_1(const std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, false, _geneticFilename, 1, 1) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x1, 0x0, 0x1};
  char readHeader[8];

//  Rcpp::Rcout << "Contructor 1.1" << std::endl;  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
//    Rcpp::Rcout << "Header:\t" << (int)readHeader[4] << '\t'
//                << (int)readHeader[5] << '\t'
//                << (int)readHeader[6] << '\t'
//                << (int)readHeader[7] << std::endl;
    if (std::memcmp(header, readHeader, 8)) {
      m_errorMessage = "File is not a binary dosage format 1.1 file";
    } else {
      m_valid = true;
//      Rcpp::Rcout << m_errorMessage << std::endl;
    }
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
      *d = *u / (double)0x7ffe;
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
  m_infile.seekg(8 + (n - 1) * m_numSubjects * sizeof(unsigned short));
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good()) {
    m_errorMessage = "Dosage read failure";
    return 1;
  }
  AssignDosages();
  return 0;
}

// ******************   CBinaryDosageFormat1_2   **********************/

CBinaryDosageFormat1_2::CBinaryDosageFormat1_2(const std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, true, _geneticFilename, 1, 2) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x1, 0x0, 0x2};
  char readHeader[8];
  
  //  Rcpp::Rcout << "Contructor 1.1" << std::endl;  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
    //    Rcpp::Rcout << "Header:\t" << (int)readHeader[4] << '\t'
    //                << (int)readHeader[5] << '\t'
    //                << (int)readHeader[6] << '\t'
    //                << (int)readHeader[7] << std::endl;
    if (std::memcmp(header, readHeader, 8)) {
      m_errorMessage = "File is not a binary dosage format 1.1 file";
    } else {
      m_valid = true;
      //      Rcpp::Rcout << m_errorMessage << std::endl;
    }
  }
  m_readBuffer.resize(2 * m_numSubjects);
}

void CBinaryDosageFormat1_2::AssignDosages() {
  unsigned short *u1, *u2;
  double *d, *p0, *p1, *p2;
  unsigned int ui;

  u1 = &m_readBuffer[0];
  u2 = u1 + m_numSubjects;
  p1 = &m_probabilities[1][0];
  p2 = &m_probabilities[2][0];
  for (ui = 0; ui < m_numSubjects; ++ui, ++u1, ++u2, ++p1, ++p2) {
    if (*u1 == 0xffff) {
      *p1 = NA_REAL;
      *p2 = NA_REAL;
    } else {
      *p1 = *u1 / (double)0xfffe;
      *p2 = *u2 / (double)0xfffe;
    }
  }
  d = &m_dosages[0];
  p0 = &m_probabilities[0][0];
  p1 = &m_probabilities[1][0];
  p2 = &m_probabilities[2][0];
  for (ui = 0; ui < m_numSubjects; ++ui, ++d, ++p0, ++p1, ++p2) {
    if (*p1 == *p1) { // Not a number (nan) check
      *d = *p1 + *p2 + *p2;
      if (*d > 2)
        *d = 2;
      *p0 = 1 - *p1 - *p2;
      if (*p0 < 0)
        *p0 = 0;
    } else {
      *d = NA_REAL;
      *p0 = NA_REAL;
    }
  }
}

int CBinaryDosageFormat1_2::GetFirst() {
  return GetSNP(1);
}

int CBinaryDosageFormat1_2::GetNext() {
  if (!m_valid)
    return 1;
  m_infile.read((char *)&m_readBuffer[0], 2 * m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}

int CBinaryDosageFormat1_2::GetSNP(unsigned int n) {
  if (CGeneticData::GetSNP(n))
    return 1;
  m_infile.clear();
  m_infile.seekg(8 + (n - 1) * m_numSubjects * sizeof(unsigned short));
  m_infile.read((char *)&m_readBuffer[0], 2 * m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good()) {
    m_errorMessage = "Dosage read failure";
    return 1;
  }
  AssignDosages();
  return 0;
}

// ******************   CBinaryDosageFormat2_1   **********************/

CBinaryDosageFormat2_1::CBinaryDosageFormat2_1(const std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, false, _geneticFilename, 1, 1) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x2, 0x0, 0x1};
  char readHeader[8];
  
  //  Rcpp::Rcout << "Contructor 2.1" << std::endl;  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
    //    Rcpp::Rcout << "Header:\t" << (int)readHeader[4] << '\t'
    //                << (int)readHeader[5] << '\t'
    //                << (int)readHeader[6] << '\t'
    //                << (int)readHeader[7] << std::endl;
    if (std::memcmp(header, readHeader, 8)) {
      m_errorMessage = "File is not a binary dosage format 1.2 file";
    } else {
      m_valid = true;
      //      Rcpp::Rcout << m_errorMessage << std::endl;
    }
  }
  m_readBuffer.resize(m_numSubjects);
}

void CBinaryDosageFormat2_1::AssignDosages() {
  unsigned short *u;
  double *d;
  unsigned int ui;
  
  for (ui = 0, u = &m_readBuffer[0], d = &m_dosages[0]; ui < m_numSubjects; ++ui, ++u, ++d) {
    if (*u > 20000)
      *d = NA_REAL;
    else
      *d = *u / 10000.;
  }
}

int CBinaryDosageFormat2_1::GetFirst() {
  return GetSNP(1);
}

int CBinaryDosageFormat2_1::GetNext() {
  if (!m_valid)
    return 1;
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}

int CBinaryDosageFormat2_1::GetSNP(unsigned int n) {
  if (CGeneticData::GetSNP(n))
    return 1;
  m_infile.clear();
  m_infile.seekg(8 + (n - 1) * m_numSubjects * sizeof(unsigned short));
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good()) {
    m_errorMessage = "Dosage read failure";
    return 1;
  }
  AssignDosages();
  return 0;
}

// ******************   CBinaryDosageFormat2_2   **********************/

CBinaryDosageFormat2_2::CBinaryDosageFormat2_2(const std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, true, _geneticFilename, 1, 2) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x2, 0x0, 0x2};
  char readHeader[8];
  
  //  Rcpp::Rcout << "Contructor 1.1" << std::endl;  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
    //    Rcpp::Rcout << "Header:\t" << (int)readHeader[4] << '\t'
    //                << (int)readHeader[5] << '\t'
    //                << (int)readHeader[6] << '\t'
    //                << (int)readHeader[7] << std::endl;
    if (std::memcmp(header, readHeader, 8)) {
      m_errorMessage = "File is not a binary dosage format 2.1 file";
    } else {
      m_valid = true;
      //      Rcpp::Rcout << m_errorMessage << std::endl;
    }
  }
  m_readBuffer.resize(2 * m_numSubjects);
}

void CBinaryDosageFormat2_2::AssignDosages() {
  unsigned short *u1, *u2;
  double *d, *p0, *p1, *p2;
  unsigned int ui;
  
  u1 = &m_readBuffer[0];
  u2 = u1 + m_numSubjects;
  p1 = &m_probabilities[1][0];
  p2 = &m_probabilities[2][0];
  for (ui = 0; ui < m_numSubjects; ++ui, ++u1, ++u2, ++p1, ++p2) {
    if (*u1 > 20000) {
      *p1 = NA_REAL;
      *p2 = NA_REAL;
    } else {
      *p1 = *u1 / 10000.;
      *p2 = *u2 / 10000.;
    }
  }
  d = &m_dosages[0];
  p0 = &m_probabilities[0][0];
  p1 = &m_probabilities[1][0];
  p2 = &m_probabilities[2][0];
  for (ui = 0; ui < m_numSubjects; ++ui, ++d, ++p0, ++p1, ++p2) {
    if (*p1 == *p1) { // Not a number (nan) check
      *d = *p1 + *p2 + *p2;
      if (*d > 2)
        *d = 2;
      *p0 = 1 - *p1 - *p2;
      if (*p0 < 0)
        *p0 = 0;
    } else {
      *d = NA_REAL;
      *p0 = NA_REAL;
    }
  }
}

int CBinaryDosageFormat2_2::GetFirst() {
  return GetSNP(1);
}

int CBinaryDosageFormat2_2::GetNext() {
  if (!m_valid)
    return 1;
  m_infile.read((char *)&m_readBuffer[0], 2 * m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}

int CBinaryDosageFormat2_2::GetSNP(unsigned int n) {
  if (CGeneticData::GetSNP(n))
    return 1;
  m_infile.clear();
  m_infile.seekg(8 + (n - 1) * m_numSubjects * sizeof(unsigned short));
  m_infile.read((char *)&m_readBuffer[0], 2 * m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good()) {
    m_errorMessage = "Dosage read failure";
    return 1;
  }
  AssignDosages();
  return 0;
}

// ******************   CBinaryDosageFormat3_1   **********************/

CBinaryDosageFormat3_1::CBinaryDosageFormat3_1(const std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, false, _geneticFilename, 1, 1) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x3, 0x0, 0x1};
  char readHeader[8];
  
  //  Rcpp::Rcout << "Contructor 2.1" << std::endl;  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
    //    Rcpp::Rcout << "Header:\t" << (int)readHeader[4] << '\t'
    //                << (int)readHeader[5] << '\t'
    //                << (int)readHeader[6] << '\t'
    //                << (int)readHeader[7] << std::endl;
    if (std::memcmp(header, readHeader, 8)) {
      m_errorMessage = "File is not a binary dosage format 2.2 file";
    } else {
      m_valid = true;
      //      Rcpp::Rcout << m_errorMessage << std::endl;
    }
  }
  m_readBuffer.resize(m_numSubjects);
}

void CBinaryDosageFormat3_1::AssignDosages() {
  unsigned short *u;
  double *d;
  unsigned int ui;
  
  for (ui = 0, u = &m_readBuffer[0], d = &m_dosages[0]; ui < m_numSubjects; ++ui, ++u, ++d) {
    if (*u > 20000)
      *d = NA_REAL;
    else
      *d = *u / 10000.;
  }
}

int CBinaryDosageFormat3_1::GetFirst() {
  return GetSNP(1);
}

int CBinaryDosageFormat3_1::GetNext() {
  if (!m_valid)
    return 1;
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}

int CBinaryDosageFormat3_1::GetSNP(unsigned int n) {
  if (CGeneticData::GetSNP(n))
    return 1;
  m_infile.clear();
  m_infile.seekg(12 + (n - 1) * m_numSubjects * sizeof(unsigned short));
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good()) {
    m_errorMessage = "Dosage read failure";
    return 1;
  }
  AssignDosages();
  return 0;
}

// ******************   CBinaryDosageFormat3_2   **********************/

CBinaryDosageFormat3_2::CBinaryDosageFormat3_2(const std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, true, _geneticFilename, 1, 2) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x3, 0x0, 0x2};
  char readHeader[8];
  
  //  Rcpp::Rcout << "Contructor 1.1" << std::endl;  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
    //    Rcpp::Rcout << "Header:\t" << (int)readHeader[4] << '\t'
    //                << (int)readHeader[5] << '\t'
    //                << (int)readHeader[6] << '\t'
    //                << (int)readHeader[7] << std::endl;
    if (std::memcmp(header, readHeader, 8)) {
      m_errorMessage = "File is not a binary dosage format 3.2 file";
    } else {
      m_valid = true;
      //      Rcpp::Rcout << m_errorMessage << std::endl;
    }
  }
  m_readBuffer.resize(4 * m_numSubjects);
}

void CBinaryDosageFormat3_2::AssignDosages() {
  unsigned short *u1, *u2;
  double *d, *p0, *p1, *p2;
  unsigned int ui;
  
  u1 = &m_readBuffer[0];
  u2 = u1 + m_numSubjects;
  d = &m_dosages[0];
  p0 = &m_probabilities[0][0];
  p1 = &m_probabilities[1][0];
  p2 = &m_probabilities[2][0];
  for (ui = 0; ui < m_numSubjects; ++ui, ++u1, ++d, ++p0, ++p1, ++p2) {
    *d = (*u1 & 0x7fff) / 10000.;
    if (*d > 2.) {
      *d = NA_REAL;
      *p0 = NA_REAL;
      *p1 = NA_REAL;
      *p2 = NA_REAL;
    } else if (*u1 & 0x8000) {
      *p1 = (*u2 & 0x7fff) / 10000.;
      if (*u2 & 0x8000) {
        ++u2;
        *p0 = *u2 / 10000.;
        ++u2;
        *p2 = *u2 / 10000.;
      } else {
        *p2 = (*d - *p1) / 2;
        *p0 = 1 - *p2 - *p1;
      }
      ++u2;
    } else {
      if (*d > 1.) {
        *p0 = 0;
        *p2 = *d - 1.;
        *p1 = 1. - *p2;
      } else {
        *p2 = 0;
        *p1 = *d;
        *p0 = 1 - *p1;
      }
    }
  }
}

int CBinaryDosageFormat3_2::GetFirst() {
  return GetSNP(1);
}

int CBinaryDosageFormat3_2::GetNext() {
  int snpSize;

  if (!m_valid)
    return 1;
  m_infile.read((char *)&snpSize, sizeof(int));
  m_infile.read((char *)&m_readBuffer[0], snpSize);
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}

int CBinaryDosageFormat3_2::GetSNP(unsigned int n) {
  int snpSize;
  unsigned int ui;
  
  if (CGeneticData::GetSNP(n))
    return 1;
  m_infile.clear();
  m_infile.seekg(12);
  for (ui = 1; ui < n; ++ui) {
    m_infile.read((char *)&snpSize, sizeof(int));
    m_infile.seekg(snpSize, std::ios_base::cur);
  }
  m_infile.read((char *)&snpSize, sizeof(int));
  m_infile.read((char *)&m_readBuffer[0], snpSize);
  if (!m_infile.good()) {
    m_errorMessage = "Dosage read failure";
    return 1;
  }
  AssignDosages();
  return 0;
}

// ******************   CBinaryDosageFormat4_1   **********************/

CBinaryDosageFormat4_1::CBinaryDosageFormat4_1(const std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, true, _geneticFilename, 1, 1) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x4, 0x0, 0x1};
  char readHeader[8];
  unsigned int beginDosages;
  
  m_firstSNPpos = 0;
  //  Rcpp::Rcout << "Contructor 1.1" << std::endl;  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
    //    Rcpp::Rcout << "Header:\t" << (int)readHeader[4] << '\t'
    //                << (int)readHeader[5] << '\t'
    //                << (int)readHeader[6] << '\t'
    //                << (int)readHeader[7] << std::endl;
    if (std::memcmp(header, readHeader, 8)) {
      m_errorMessage = "File is not a binary dosage format 4.2 file";
    } else {
      m_infile.seekg(36);
      m_infile.read((char *)&beginDosages, sizeof(unsigned int));
      if (!m_infile.good()) {
        m_errorMessage = "Error finding start of dosage data";
      } else {
        m_firstSNPpos = beginDosages;
        //        Rcpp::Rcout << "First dosage postion:\t" << m_firstSNPpos << std::endl;
        m_valid = true;
      }
      //      Rcpp::Rcout << m_errorMessage << std::endl;
    }
  }
  m_readBuffer.resize(m_numSubjects);
}

void CBinaryDosageFormat4_1::AssignDosages() {
  unsigned short *u;
  double *d;
  unsigned int ui;
  
  for (ui = 0, u = &m_readBuffer[0], d = &m_dosages[0]; ui < m_numSubjects; ++ui, ++u, ++d) {
    if (*u > 20000)
      *d = NA_REAL;
    else
      *d = *u / 10000.;
  }
}

int CBinaryDosageFormat4_1::GetFirst() {
  return GetSNP(1);
}

int CBinaryDosageFormat4_1::GetNext() {
  if (!m_valid)
    return 1;
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}

int CBinaryDosageFormat4_1::GetSNP(unsigned int n) {
  if (CGeneticData::GetSNP(n))
    return 1;
  m_infile.clear();
  m_infile.seekg(m_firstSNPpos + (n - 1) * m_numSubjects * sizeof(unsigned short));
  m_infile.read((char *)&m_readBuffer[0], m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good()) {
    m_errorMessage = "Dosage read failure";
    return 1;
  }
  AssignDosages();
  return 0;
}

// ******************   CBinaryDosageFormat4_2   **********************/

CBinaryDosageFormat4_2::CBinaryDosageFormat4_2(const std::string &_geneticFilename, const int _numSubjects, const int _numSNPs)
  : CBinaryDosage(_numSubjects, _numSNPs, true, _geneticFilename, 1, 2) {
  const char header[8] = {'b', 'o', 's', 'e', 0x0, 0x4, 0x0, 0x2};
  char readHeader[8];
  unsigned int beginDosages;

  m_firstSNPpos = 0;
  //  Rcpp::Rcout << "Contructor 1.1" << std::endl;  
  if (m_infile.good()) {
    m_infile.read(readHeader, 8);
    //    Rcpp::Rcout << "Header:\t" << (int)readHeader[4] << '\t'
    //                << (int)readHeader[5] << '\t'
    //                << (int)readHeader[6] << '\t'
    //                << (int)readHeader[7] << std::endl;
    if (std::memcmp(header, readHeader, 8)) {
      m_errorMessage = "File is not a binary dosage format 4.2 file";
    } else {
      m_infile.seekg(36);
      m_infile.read((char *)&beginDosages, sizeof(unsigned int));
      if (!m_infile.good()) {
        m_errorMessage = "Error finding start of dosage data";
      } else {
        m_firstSNPpos = beginDosages;
//        Rcpp::Rcout << "First dosage postion:\t" << m_firstSNPpos << std::endl;
        m_valid = true;
      }
      //      Rcpp::Rcout << m_errorMessage << std::endl;
    }
  }
  m_readBuffer.resize(4 * m_numSubjects);
}

void CBinaryDosageFormat4_2::AssignDosages() {
  unsigned short *u1, *u2;
  double *d, *p0, *p1, *p2;
  unsigned int ui;
  
  u1 = &m_readBuffer[0];
  u2 = u1 + m_numSubjects;
  d = &m_dosages[0];
  p0 = &m_probabilities[0][0];
  p1 = &m_probabilities[1][0];
  p2 = &m_probabilities[2][0];
  for (ui = 0; ui < m_numSubjects; ++ui, ++u1, ++d, ++p0, ++p1, ++p2) {
    *d = (*u1 & 0x7fff) / 10000.;
    if (*d > 2.) {
      *d = NA_REAL;
      *p0 = NA_REAL;
      *p1 = NA_REAL;
      *p2 = NA_REAL;
    } else if (*u1 & 0x8000) {
      *p1 = (*u2 & 0x7fff) / 10000.;
      if (*u2 & 0x8000) {
        ++u2;
        *p0 = *u2 / 10000.;
        ++u2;
        *p2 = *u2 / 10000.;
      } else {
        *p2 = (*d - *p1) / 2;
        *p0 = 1 - *p2 - *p1;
      }
      ++u2;
    } else {
      if (*d > 1.) {
        *p0 = 0;
        *p2 = *d - 1.;
        *p1 = 1. - *p2;
      } else {
        *p2 = 0;
        *p1 = *d;
        *p0 = 1 - *p1;
      }
    }
  }
}

int CBinaryDosageFormat4_2::GetFirst() {
  return GetSNP(1);
}

int CBinaryDosageFormat4_2::GetNext() {
  int snpSize;
  
  if (!m_valid)
    return 1;
  m_infile.read((char *)&snpSize, sizeof(int));
  m_infile.read((char *)&m_readBuffer[0], snpSize);
  if (!m_infile.good())
    return 1;
  AssignDosages();
  return 0;
}

int CBinaryDosageFormat4_2::GetSNP(unsigned int n) {
  int snpSize;
  unsigned int ui;
  
  if (CGeneticData::GetSNP(n))
    return 1;
  m_infile.clear();
  m_infile.seekg(m_firstSNPpos);
  for (ui = 1; ui < n; ++ui) {
    m_infile.read((char *)&snpSize, sizeof(int));
    m_infile.seekg(snpSize, std::ios_base::cur);
  }
  m_infile.read((char *)&snpSize, sizeof(int));
  m_infile.read((char *)&m_readBuffer[0], snpSize);
  if (!m_infile.good()) {
    m_errorMessage = "Dosage read failure";
    return 1;
  }
  AssignDosages();
  return 0;
}
