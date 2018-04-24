#ifndef IMP2TOBD
#define IMP2TOBD 1

#ifndef IMPUTE2_H
#include "Impute2.h"
#endif

int Imp2toBDC(const Rcpp::List &imp2Info, const std::string &filename, int format, int subformat);

class CWriteBD {
protected:
  std::ofstream m_outfile;
  int m_numSubjects;
  int m_numSNPs;
  std::vector<unsigned short> m_d, m_p0, m_p1, m_p2;
  std::vector<unsigned short> m_additional;
  
  CWriteBD(const std::string &_filename, int _numSub, int _numSNPs);
public:
  virtual ~CWriteBD() { m_outfile.close(); }

  virtual int WriteHeader();
  virtual int WriteSNP(const double *_d);
  virtual int WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2);
};

class CWriteBD11 : public CWriteBD {
public:
  CWriteBD11(const std::string &_filename, int _numSub, int _numSNPs);
  
  virtual int WriteHeader();
  virtual int WriteSNP(const double *_d);
  virtual int WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2);
};

class CWriteBD12 : public CWriteBD {
public:
  CWriteBD12(const std::string &_filename, int _numSub, int _numSNPs);
  
  virtual int WriteHeader();
  virtual int WriteSNP(const double *_d);
  virtual int WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2);
};

class CWriteBD21 : public CWriteBD {
public:
  CWriteBD21(const std::string &_filename, int _numSub, int _numSNPs);
  
  virtual int WriteHeader();
  virtual int WriteSNP(const double *_d);
  virtual int WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2);
};

class CWriteBD22 : public CWriteBD {
public:
  CWriteBD22(const std::string &_filename, int _numSub, int _numSNPs);
  
  virtual int WriteHeader();
  virtual int WriteSNP(const double *_d);
  virtual int WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2);
};

class CWriteBD31 : public CWriteBD {
public:
  CWriteBD31(const std::string &_filename, int _numSub, int _numSNPs);
  
  virtual int WriteHeader();
  virtual int WriteSNP(const double *_d);
  virtual int WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2);
};

class CWriteBD32 : public CWriteBD {
public:
  CWriteBD32(const std::string &_filename, int _numSub, int _numSNPs);
  
  virtual int WriteHeader();
  virtual int WriteSNP(const double *_d);
  virtual int WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2);
};

class CWriteBD41 : public CWriteBD {
public:
  CWriteBD41(const std::string &_filename, int _numSub, int _numSNPs);
  
  virtual int WriteHeader();
  virtual int WriteSubjectAndSNPInfo(const Rcpp::List &imp2Info);
  virtual int WriteSNP(const double *_d);
  virtual int WriteSNP(const double *_d, const double *_p0, const double *_p1, const double *_p2);
};

#endif
