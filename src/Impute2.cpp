#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <Rcpp.h>
#include "Impute2.h"

CGeneticData *OpenImpute2File(const Rcpp::List &geneticInfo) {
  CGeneticData *geneticData;
  std::string filename;
  int numSubjects, numSNPs;
  bool header;
  std::vector<int> snpCol;
  int startCol;
  int format;
  char sep = '\t';
  int i;
  
  filename = Rcpp::as<std::string>(geneticInfo["filename"]);
  numSubjects = (int)geneticInfo["numSubjects"];
  numSNPs = (int)geneticInfo["numSNPs"];
  header = (bool)geneticInfo["header"];
  snpCol = Rcpp::as<std::vector<int> >(geneticInfo["snpCol"]);
  startCol = (int)geneticInfo["startCol"];
  format = (int)geneticInfo["format"];
  geneticData = new CImpute2(numSubjects, numSNPs, FALSE, (format != 1),
                             filename, header, snpCol, startCol, format , sep);
  return geneticData;
}

CImpute2::CImpute2(int _numSubjects, int _numSNPs, bool _measured, bool _geneticProbabilities,
                   const std::string &_filename, bool _header, const std::vector<int> &_snpCol, int _startCol,
                   int _format, char _sep) : CGeneticData(_numSubjects, _numSNPs, _measured, _geneticProbabilities) {
  m_filename = _filename;
  m_header = _header;
  m_snpCol = _snpCol;
  m_startCol = _startCol;
  m_format = _format;
  if (m_doseData)
    delete [] m_doseData;
  m_doseData = new double[m_numSubjects];
  m_dosages = m_doseData;
  if (m_format != 1) {
    m_probabilities[0] = m_dosages + m_numSubjects;
    m_probabilities[1] = m_probabilities[0] + m_numSubjects;
    m_probabilities[2] = m_probabilities[1] + m_numSubjects;
  } else {
    m_probabilities[0] = NULL;
    m_probabilities[1] = NULL;
    m_probabilities[2] = NULL;
  }
  m_sep = _sep;
  
  m_infile.open(m_filename);
  if (!m_infile.good()) {
    m_errorMessage = "Unable to open file";
    m_valid = false;
    return;
  }
  m_valid = true;
}

int CImpute2::GetFirst() {
  std::string junk;
  
  m_infile.clear();
  m_infile.seekg(0);
  if (m_header)
    std::getline(m_infile, junk);
  return GetNext();
}

int CImpute2::GetNext() {
  std::string readline;
  std::string junk;
  std::istringstream istring;
  int i, j;
  
  std::getline(m_infile, readline);
  if (!m_infile.good()) {
    m_errorMessage = "Error reading line from file";
    return 1;
  }
  istring.str(readline);
  for (i = 1; i < m_startCol; ++i)
    istring >> junk;
  if (!istring.good())
    m_errorMessage = "Error reading line from file";
  if (m_format == 1) {
    for (i = 0; i < m_numSubjects; ++i)
      istring >> m_dosages[i];
  } else {
    for (i = 0; i < m_numSubjects; ++i) {
      for (j = 0; j < m_format; ++j)
        istring >> m_probabilities[j][i];
      if (m_format == 2) {
        m_probabilities[2][i] = 1. - m_probabilities[0][i] - m_probabilities[1][i];
        if (m_probabilities[2][i] < 0.)
          m_probabilities[2][i] = 0.;
      }
      m_dosages[i] = m_probabilities[1][i] + m_probabilities[2][i] + m_probabilities[2][i];
      if (m_dosages[i] > 2.)
        m_dosages[i] = 2.;
    }
  }
  if (m_infile.bad()) {
    m_errorMessage = "Error reading line from file";
    return 1;
  }
  
  return 0;
}

int CImpute2::GetSNP(unsigned int n) {
  std::string junk;
  int ui;
  
  if (CGeneticData::GetSNP(n))
    return 1;
  
  if (m_header)
    std::getline(m_infile, junk);
  for (ui = 0; ui < n; ++ui)
    std::getline(m_infile, junk);

  if (!m_infile.good()) {
    m_errorMessage = "Error reading from file";
    return 1;
  }
  return GetNext();
}