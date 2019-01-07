#include <iostream>
#include <fstream>
#include "Rcpp.h"

// [[Rcpp::export]]
int GetLocations(Rcpp::IntegerVector x, Rcpp::NumericVector y, double fileSize, int bufferSize) {
  int *xi;
  long long *loc;
  long long currentSectionStart;
  long long bufSize = bufferSize;
  int numSections;

  xi = &x[0];
  loc = (long long *)&y[0];
  *loc = *xi;
  ++xi;
  ++loc;
  for (int i = 1; i < x.size(); ++i, ++xi, ++loc) {
    *loc = *(loc - 1) + *xi;
  }
  *loc = fileSize;
  Rcpp::Rcout << "File size\t" << *loc << std::endl;
  loc = (long long *)&y[0];
  currentSectionStart = *loc;
  ++loc;
  numSections = 1;
  for (int i = 1; i < y.size(); ++i, ++loc) {
    if ((*loc - currentSectionStart) > bufSize) {
      ++numSections;
      currentSectionStart = *(loc - 1);
    }
  }
  
  return numSections;
}

// [[Rcpp::export]]
int WriteLocations(Rcpp::NumericVector y) {
  long long *loc;
  
  loc = (long long *)&y[0];
  for (int i = 0; i < 5; ++i, ++loc)
    Rcpp::Rcout << *loc << '\t';
  loc = (long long *)&y[0];
  Rcpp::Rcout << std::endl
              << *(loc + y.size() - 5) << '\t' << *(loc + y.size() - 4) << '\t'
              << *(loc + y.size() - 3) << '\t' << *(loc + y.size() - 2) << '\t' <<  *(loc + y.size() - 1) << std::endl;
  return 0;
}

// [[Rcpp::export]]
int GetSections(Rcpp::NumericVector y, Rcpp::IntegerVector snpSection, Rcpp::NumericVector fileLocation, Rcpp::IntegerVector snpLocation, int bufferSize) {
  long long *loc, *fileLoc;
  long long currentSectionStart;
  long long bufSize = bufferSize;
  int *snpSec, *snpLoc;
  int currentSection;
  
  loc = (long long *)&y[0];
  currentSectionStart = *loc;
  fileLoc = (long long *)&fileLocation[0];
  currentSection = 0;
  snpSec = &snpSection[0];
  snpLoc = &snpLocation[0];
  for (int i = 0; i < snpSection.size(); ++i, ++snpSec, ++snpLoc) {
    *snpLoc = (int)(*loc - currentSectionStart);
    ++loc;
    if ((*loc - currentSectionStart) > bufSize) {
      ++currentSection;
      *fileLoc = currentSectionStart;
      ++fileLoc;
      currentSectionStart = *(loc - 1);
      *snpLoc = 0;
    }
    *snpSec = currentSection;
  }
  *fileLoc = currentSectionStart;
  ++fileLoc;
  *fileLoc = *loc;
  
  return 0;
}

// [[Rcpp::export]]
int ReadSNP(Rcpp::IntegerVector &snpNumber, std::string &filename, Rcpp::IntegerVector &format,
            int numSub, int numSNPs, int bufferSize, Rcpp::IntegerVector &buffer,
            int sections, Rcpp::IntegerVector &snpSection, Rcpp::NumericVector &fileLocation,
            Rcpp::IntegerVector &snpLocation, int currentSection, Rcpp::NumericVector &values) {
  
  std::ifstream infile;
  long long *fLoc;
  
  fLoc = (long long *)&fileLocation[0];
  for (int i = 0; i < snpNumber.size(); ++i) {
    if (snpNumber[i] > numSNPs)
      return 1;
    if (snpSection[snpNumber[i]] != currentSection) {
      if (!infile.is_open()) {
        infile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
        if (!infile.good())
          return 1;
      }
    }
    infile.seekg(fLoc[snpSection[snpNumber[i]]]);
    Rcpp::Rcout << fLoc[snpSection[snpNumber[i]] + 1] - fLoc[snpSection[snpNumber[i]]] << std::endl;
    infile.read((char *)&buffer[0], fLoc[snpSection[snpNumber[i]] + 1] - fLoc[snpSection[snpNumber[i]]]);
  }
  
  return 0;
}