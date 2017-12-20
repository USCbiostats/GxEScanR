#ifndef PLINKBINARY_H
#define PLINKBINARY_H 1

#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

Rcpp::List PlinkBinaryInfo(std::string &geneticFile, std::string &mapFile, std::string &familyFile);

class CPlinkBinary : public CGeneticData {
protected:
  std::string m_geneticFile;
  std::string m_mapFile;
  std::string m_familyFile;
  
  int CheckFileSize();
public:
  CPlinkBinary(Rcpp::List &_geneticData);
  CPlinkBinary(std::string &geneticFile, std::string &mapFile, std::string &familyFile);
  virtual ~CPlinkBinary() {}
  
  virtual int CheckValidity();
  virtual int GetFirst() { return 1; }
  virtual int GetNext() { return 1; }
};

#endif
