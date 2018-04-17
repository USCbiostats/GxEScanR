#ifndef IMPUTE2_H
#define IMPUTE2_H 1
#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

class CImpute2 : public CGeneticData {
protected:
  std::ifstream m_infile;
  std::string m_filename;
  bool m_header;
  std::vector<int> m_snpCol;
  int m_startCol;
  int m_format;
  char m_sep;
public:
  CImpute2(int _numSubjects, int _numSNPs, bool _measured, bool _geneticProbabilities,
           const std::string &_filename, bool _header, const std::vector<int> &_snpCol, int _startCol,
           int _format, char _sep);
  virtual ~CImpute2() { m_infile.close(); }

  virtual int GetFirst();
  virtual int GetNext();
  virtual int GetSNP(unsigned int n);
};

#endif