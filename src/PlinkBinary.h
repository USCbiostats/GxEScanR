#ifndef PLINKBINARY_H
#define PLINKBINARY_H 1

#ifndef GENETICDATA
#include "GeneticData.h"
#endif

Rcpp::List PlinkBinaryInfo(std::string &geneticFile, std::string &mapFile, std::string &familyFile);

class CPlinkBinary : public CGeneticData {
protected:
  std::string m_geneticFile;
  std::string m_mapFile;
  std::string m_familyFile;
public:
  CPlinkBinary(std::string &geneticFile, std::string &mapFile, std::string &familyFile);
  virtual ~CPlinkBinary() {}
  virtual int CheckValidity();
};

#endif
