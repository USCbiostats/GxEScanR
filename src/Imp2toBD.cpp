#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <Rcpp.h>
#include "Imp2toBD.h"

unsigned short DtoU(double d) {
  short u1, u2;
  double d1, d2;
  
  u1 = (short)(d * 10000);
  d1 = ((double)u1) / 10000;
  if (d1 > d)
    u2 = u1 - 1;
  else
    u2 = u1 + 1;
  d2 = ((double)u2) / 100000;
  if (fabs(d2 - d) < fabs(d1 - d))
    return u2;
  return u1;
}
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
  } else if (format ==1 && subformat == 2) {
    bdFile = new CWriteBD12(filename, (int)imp2File->NumSubjects(), (int)imp2File->NumSNPs());
  } else if (format == 2 && subformat == 1) {
    bdFile = new CWriteBD21(filename, (int)imp2File->NumSubjects(), (int)imp2File->NumSNPs());
  } else if (format == 2 && subformat == 2) {
    bdFile = new CWriteBD22(filename, (int)imp2File->NumSubjects(), (int)imp2File->NumSNPs());
  } else if (format == 3 && subformat == 1) {
    bdFile = new CWriteBD31(filename, (int)imp2File->NumSubjects(), (int)imp2File->NumSNPs());
  } else if (format == 3 && subformat == 2) {
    bdFile = new CWriteBD32(filename, (int)imp2File->NumSubjects(), (int)imp2File->NumSNPs());
  } else {
    Rcpp::Rcout << "Unknown bdosage format";
    return 1;
  }
  bdFile->WriteHeader();
  Rcpp::Rcout << "Wrote header" << std::endl;

  imp2File->GetFirst();
  if (subformat == 1)
    bdFile->WriteSNP(imp2File->Dosages().data());
  else
    bdFile->WriteSNP(imp2File->Dosages().data(), imp2File->Probabilities()[0].data(), imp2File->Probabilities()[1].data(), imp2File->Probabilities()[2].data());

  for (i = 1; i < imp2File->NumSNPs(); ++i) {
    imp2File->GetNext();
    if (subformat == 1)
      bdFile->WriteSNP(imp2File->Dosages().data());
    else
      bdFile->WriteSNP(imp2File->Dosages().data(), imp2File->Probabilities()[0].data(), imp2File->Probabilities()[1].data(), imp2File->Probabilities()[2].data());
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
  m_additional.resize(4 * m_numSubjects);
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

// *******************   Format 1.1  ******************************************

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
  for (i = 0; i < m_numSubjects; ++i, ++d, ++u) {
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

// *******************   Format 1.2  ******************************************

CWriteBD12::CWriteBD12(const std::string &_filename, int _numSub, int _numSNPs) : CWriteBD(_filename, _numSub, _numSNPs) {}

int CWriteBD12::WriteHeader() {
  const char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x1, 0x0, 0x2 };
  if (CWriteBD::WriteHeader())
    return 1;
  
  m_outfile.write(header, 8);
  return 0;
}

int CWriteBD12::WriteSNP(const double *_d) {
  return 1;
}

int CWriteBD12::WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2) {
  const double *p1, *p2;
  double d1, d2;
  unsigned short *u1, *u2;
  unsigned short up1, up2;
  int i;
  
  if (CWriteBD::WriteSNP(_d, _p0, _p1, _p2))
    return 1;
  
  p1 = _p1;
  p2 = _p2;
  u1 = m_p1.data();
  u2 = m_p2.data();
  
  for (i = 0; i < m_numSubjects; ++i, ++p1, ++p2, ++u1, ++u2) {
    up1 = (short)(*p1 * 0xfffe);
    d1 = ((double)up1) / 0xfffe;
    if (d1 > *p1)
      up2 = up1 - 1;
    else
      up2 = up1 + 1;
    d2 = ((double)up2) / 0xfffe;
    if (fabs(d2 - *p1) < fabs(d1 - *p1))
      *u1 = up2;
    else
      *u1 = up1;

    up1 = (short)(*p2 * 0xfffe);
    d1 = ((double)up2) / 0xfffe;
    if (d1 > *p2)
      up2 = up1 - 1;
    else
      up2 = up1 + 1;
    d2 = ((double)up2) / 0xfffe;
    if (fabs(d2 - *p2) < fabs(d1 - *p2))
      *u2 = up2;
    else
      *u2 = up1;
    
    //    if (i < 10)
    //      Rcpp::Rcout << *d << '\t' << *u << '\t' << u1 << '\t' << d1 << '\t' << u2 << '\t' << d2 << std::endl;
  }
  u1 = m_p1.data();
  u2 = m_p2.data();
  m_outfile.write((char *)u1, sizeof(unsigned short) * m_numSubjects);
  m_outfile.write((char *)u2, sizeof(unsigned short) * m_numSubjects);
  
  return 0;
}

// *******************   Format 2.1  ******************************************

CWriteBD21::CWriteBD21(const std::string &_filename, int _numSub, int _numSNPs) : CWriteBD(_filename, _numSub, _numSNPs) {}

int CWriteBD21::WriteHeader() {
  const char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x2, 0x0, 0x1 };
  if (CWriteBD::WriteHeader())
    return 1;
  
  m_outfile.write(header, 8);
  return 0;
}

int CWriteBD21::WriteSNP(const double *_d) {
  const double *d;
  double d1, d2;
  unsigned short *u;
  unsigned short u1, u2;
  int i;
  
  if (CWriteBD::WriteSNP(_d))
    return 1;
  
  d = _d;
  u = m_d.data();
  for (i = 0; i < m_numSubjects; ++i, ++d, ++u) {
    u1 = (short)(*d * 10000);
    d1 = ((double)u1) / 10000;
    if (d1 > *d)
      u2 = u1 - 1;
    else
      u2 = u1 + 1;
    d2 = ((double)u2) / 100000;
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

int CWriteBD21::WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2) {
  if (CWriteBD::WriteSNP(_d, _p0, _p1, _p2))
    return 1;
  return 0;
}

// *******************   Format 2.2  ******************************************

CWriteBD22::CWriteBD22(const std::string &_filename, int _numSub, int _numSNPs) : CWriteBD(_filename, _numSub, _numSNPs) {}

int CWriteBD22::WriteHeader() {
  const char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x2, 0x0, 0x2 };
  if (CWriteBD::WriteHeader())
    return 1;
  
  m_outfile.write(header, 8);
  return 0;
}

int CWriteBD22::WriteSNP(const double *_d) {
  return 1;
}

int CWriteBD22::WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2) {
  const double *p1, *p2;
  double d1, d2;
  unsigned short *u1, *u2;
  unsigned short up1, up2;
  int i;
  
  if (CWriteBD::WriteSNP(_d, _p0, _p1, _p2))
    return 1;
  
  p1 = _p1;
  p2 = _p2;
  u1 = m_p1.data();
  u2 = m_p2.data();
  
  for (i = 0; i < m_numSubjects; ++i, ++p1, ++p2, ++u1, ++u2) {
    up1 = (short)(*p1 * 10000);
    d1 = ((double)up1) / 10000;
    if (d1 > *p1)
      up2 = up1 - 1;
    else
      up2 = up1 + 1;
    d2 = ((double)up2) / 10000;
    if (fabs(d2 - *p1) < fabs(d1 - *p1))
      *u1 = up2;
    else
      *u1 = up1;
    
    up1 = (short)(*p2 * 10000);
    d1 = ((double)up2) / 10000;
    if (d1 > *p2)
      up2 = up1 - 1;
    else
      up2 = up1 + 1;
    d2 = ((double)up2) / 10000;
    if (fabs(d2 - *p2) < fabs(d1 - *p2))
      *u2 = up2;
    else
      *u2 = up1;
    
    //    if (i < 10)
    //      Rcpp::Rcout << *d << '\t' << *u << '\t' << u1 << '\t' << d1 << '\t' << u2 << '\t' << d2 << std::endl;
  }
  u1 = m_p1.data();
  u2 = m_p2.data();
  m_outfile.write((char *)u1, sizeof(unsigned short) * m_numSubjects);
  m_outfile.write((char *)u2, sizeof(unsigned short) * m_numSubjects);
  
  return 0;
}

// *******************   Format 3.1  ******************************************

CWriteBD31::CWriteBD31(const std::string &_filename, int _numSub, int _numSNPs) : CWriteBD(_filename, _numSub, _numSNPs) {}

int CWriteBD31::WriteHeader() {
  const char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x3, 0x0, 0x1 };
  if (CWriteBD::WriteHeader())
    return 1;
  
  m_outfile.write(header, 8);
  
  m_outfile.write((char *)&m_numSubjects, sizeof(int));
  return 0;
}

int CWriteBD31::WriteSNP(const double *_d) {
  const double *d;
  unsigned short *u;
  int i;
  
  if (CWriteBD::WriteSNP(_d))
    return 1;
  
  d = _d;
  u = m_d.data();
  for (i = 0; i < m_numSubjects; ++i, ++d, ++u) {
    *u = DtoU(*d);
  }
  u = m_d.data();
  m_outfile.write((char *)u, sizeof(unsigned short) * m_numSubjects);
  
  return 0;
}

int CWriteBD31::WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2) {
  if (CWriteBD::WriteSNP(_d, _p0, _p1, _p2))
    return 1;
  return 0;
}

// *******************   Format 3.2  ******************************************

CWriteBD32::CWriteBD32(const std::string &_filename, int _numSub, int _numSNPs) : CWriteBD(_filename, _numSub, _numSNPs) {}

int CWriteBD32::WriteHeader() {
  const char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x3, 0x0, 0x2 };
  if (CWriteBD::WriteHeader())
    return 1;
  
  m_outfile.write(header, 8);
  
  m_outfile.write((char *)&m_numSubjects, sizeof(int));
  return 0;
}

int CWriteBD32::WriteSNP(const double *_d) {
  if (CWriteBD::WriteSNP(_d))
    return 1;

  return 0;
}

int CWriteBD32::WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2) {
  unsigned short *u1, *u2;
  const double *d, *p0, *p1, *p2;
  int snpSize;
  int i;
  
  if (CWriteBD::WriteSNP(_d, _p0, _p1, _p2))
    return 1;
  
  snpSize = m_numSubjects;
  d = _d;
  p0 = _p0;
  p1 = _p1;
  p2 = _p2;
  u1 = m_additional.data();
  u2 = u1 + m_numSubjects;
  for (i = 0; i < m_numSubjects; ++i, ++d, ++p0, ++p1, ++p2, ++u1) {
    *u1 = DtoU(*d);
//    if (i < 4)
//      Rcpp::Rcout << fabs(*d - (*p1 + *p2 + *p2)) << '\t' << fabs(1 - *p0 - *p1 - *p2) << std::endl;
    if ((fabs(*d - (*p1 + *p2 + *p2)) > 1e-10) || (fabs(1 - *p0 - *p1 - *p2) > 1e-10)) {
//      if (i < 4)
//        Rcpp::Rcout << *d << '\t' << *p0 << '\t' << *p1 << '\t' << *p2 << std::endl;
      *u1 |= 0x8000;
      *u2 = DtoU(*p1);
      *u2 |= 0x8000;
      ++u2;
      *u2 = DtoU(*p0);
      ++u2;
      *u2 = DtoU(*p2);
      ++u2;
      snpSize += 3;
    } else if (*p1 != 0 && *p2 != 0) {
      *u1 |= 0x8000;
      *u2 = DtoU(*p1);
      ++u2;
      ++snpSize;
    }
  }
//  Rcpp::Rcout << snpSize << std::endl;
  snpSize *= sizeof(unsigned short);
  u1 = m_additional.data();
  m_outfile.write((char *)&snpSize, sizeof(int));
  m_outfile.write((char *)u1, snpSize);
  
  return 0;
}
