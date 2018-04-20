#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <Rcpp.h>
#include "Imp2toBD.h"

//' Function to convert and Impute 2 style file to a binary dosage file
//' 
//' Function to convert and Impute 2 style file to a binary dosage file
//' 
//' @param imp2Info
//' List containing information of Impute 2 file, returned from GetI2Info
//' @param filename
//' Name of binary dosage file to create
//' @param format
//' Primary format of binary dosage file
//' @param subformat
//' Secondary format of binary dosage file
//' @return
//' 0 success
//' 1 failure
//' @export
// [[Rcpp::export]]
int Imp2toBDC(const Rcpp::List &imp2Info, const std::string &filename, int format, int subformat) {
  CWriteBD *bdFile = NULL;
  CGeneticData *imp2File = NULL;
  int i;

  imp2File = OpenImpute2File(imp2Info);
  
  if (format == 1 && subformat == 1) {
    bdFile = new CWriteBD11(filename, (int)imp2File->NumSubjects(), (int)imp2File->NumSNPs());
  }
  bdFile->WriteHeader();
  imp2File->GetFirst();
  if (subformat == 1)
    bdFile->WriteSNP(imp2File->Dosages().data());
  for (i = 1; i < imp2File->NumSubjects(); ++i) {
    imp2File->GetNext();
    if (subformat == 1)
      bdFile->WriteSNP(imp2File->Dosages().data());
  }
  
  if (bdFile)
    delete bdFile;
  if (imp2File)
    delete imp2File;
  return 0;
}

CWriteBD::CWriteBD(const std::string &_filename, int _numSub, int _numSNPs) {
  Rcpp::Rcout << "Opening file" << std::endl;
  m_outfile.open(_filename.c_str(), std::ios_base::out | std::ios_base::binary);
  m_numSubjects = _numSub;
  m_numSNPs = _numSNPs;
  m_d.resize(m_numSubjects);
  m_p0.resize(m_numSubjects);
  m_p1.resize(m_numSubjects);
  m_p2.resize(m_numSubjects);
}

int CWriteBD::WriteHeader() {
  if (!m_outfile.good())
    return 1;
  return 0;
}

int CWriteBD::WriteSNP(const double *_d) {
  if (!m_outfile.good())
    return 1;
  return 0;
}

int CWriteBD::WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2) {
  if (!m_outfile.good())
    return 1;
  return 0;
}

CWriteBD11::CWriteBD11(const std::string &_filename, int _numSub, int _numSNPs) : CWriteBD(_filename, _numSub, _numSNPs) {}

int CWriteBD11::WriteHeader() {
  const char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x1, 0x0, 0x1 };
  if (CWriteBD::WriteHeader())
    return 1;
  
  m_outfile.write(header, 8);
  return 0;
}

int CWriteBD11::WriteSNP(const double *_d) {
  const double *d;
  double d1, d2;
  unsigned short *u;
  unsigned short u1, u2;
  int i;
  
  if (CWriteBD::WriteSNP(_d))
    return 1;
  
  d = _d;
  u = m_d.data();
  for (i = 0; i < m_numSubjects; ++i, ++d) {
    u1 = (short)(*d * 0x7ffe);
    d1 = ((double)u1) / 0x7ffe;
    if (d1 > *d)
      u2 = u1 - 1;
    else
      u2 = u1 + 1;
    d2 = ((double)u2) / 0x7ffe;
    if (fabs(d2 - *d) < fabs(d1 - *d))
      *u = u2;
    else
      *u = u1;
//    if (i < 10)
//      Rcpp::Rcout << *d << '\t' << *u << '\t' << u1 << '\t' << d1 << '\t' << u2 << '\t' << d2 << std::endl;
  }
  u = m_d.data();
  m_outfile.write((char *)u, sizeof(unsigned short) * m_numSubjects);
  
  return 0;
}

int CWriteBD11::WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2) {
  if (CWriteBD::WriteSNP(_d, _p0, _p1, _p2))
    return 1;
  return 0;
}
