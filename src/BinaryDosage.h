#ifndef BINARYDOSAGE_H
#define BINARYDOSAGE_H 1

#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

Rcpp::List BinaryDosageInfo(std::string &geneticFile, std::string &mapFile, std::string &familyFile);

int GetBinaryDosageFormat(std::string &binaryDosageFile, int &format, int &version, std::string &errorMessage);

#endif