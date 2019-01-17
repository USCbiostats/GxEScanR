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
  //  Rcpp::Rcout << "File size\t" << *loc << std::endl;
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


int ReadSNP42(char *bptr, int numSub, Rcpp::NumericVector &dosage, Rcpp::NumericVector &p0,
              Rcpp::NumericVector &p1, Rcpp::NumericVector &p2) {
  short *ud, *ue;
  double *pd, *pp0, *pp1, *pp2;
  
  ud = (short *)(bptr + sizeof(int));
  ue = ud + numSub;
  //  Rcpp::Rcout << "Number of bytes to read\t" << *x << std::endl;
  //  Rcpp::Rcout << "First 4 short values\t" << ud[0] << '\t' << ud[1] << '\t' << ud[2] << '\t' << ud[3] << std::endl;
  //  Rcpp::Rcout << "First extra value\t" << ue[0] << std::endl;
  pd = &dosage[0];
  pp0 = &p0[0];
  pp1 = &p1[0];
  pp2 = &p2[0];
  for (int i = 0; i < numSub; ++i, ++ud, ++pd, ++pp0, ++pp1, ++pp2) {
    if ((*ud & 0x8000) == 0x8000) {
      *pd = (*ud & 0x7fff) / 10000.;
      if ((*ue & 0x8000) == 0x8000) {
        *pp1 = (*ue & 0x7fff) / 10000.;
        ++ue;
        *pp0 = *ue / 10000.;
        ++ue;
        *pp2 = *ue / 10000.;
      } else {
        *pp1 = *ue / 10000.;
        *pp2 = (*pd - *pp1) * 0.5;
        *pp0 = 1. - *pp1 - *pp2;
      }
      ++ue;
    } else {
      *pd = *ud / 10000.;
      if (*pd > 1.) {
        *pp0 = 0.;
        *pp1 = *pd - 1.;
        *pp2 = 1. - *pp1;
      } else {
        *pp1 = *pd;
        *pp0 = 1. - *pp1;
        *pp2 = 0.;
      }
    }
  }
  return 0;
}
// [[Rcpp::export]]
int ReadSNP(Rcpp::IntegerVector &snpNumber, Rcpp::IntegerVector &subjectNumber,
            std::string &filename, Rcpp::IntegerVector &format,
            int numSub, int numSNPs, int bufferSize, Rcpp::IntegerVector &buffer,
            int sections, Rcpp::IntegerVector &snpSection, Rcpp::NumericVector &fileLocation,
            Rcpp::IntegerVector &snpLocation, int currentSection, Rcpp::NumericVector &dosage,
            Rcpp::NumericVector &p0, Rcpp::NumericVector &p1, Rcpp::NumericVector &p2,
            Rcpp::NumericMatrix &values) {
  
  std::ifstream infile;
  long long *fLoc;
  char *bptr;
  
  bptr = (char *)&buffer[0];
  fLoc = (long long *)&fileLocation[0];
  for (int i = 0; i < snpNumber.size(); ++i) {
    if (snpNumber[i] > numSNPs)
      return 1;
    if (snpSection[snpNumber[i] - 1] != currentSection) {
      if (!infile.is_open()) {
        infile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
        if (!infile.good())
          return 1;
      }
    }
    infile.seekg(fLoc[snpSection[snpNumber[i]]]);
    //    Rcpp::Rcout << "Size to read in\t" << fLoc[snpSection[snpNumber[i]] + 1] - fLoc[snpSection[snpNumber[i]]] << std::endl;
    infile.read((char *)&buffer[0], fLoc[snpSection[snpNumber[i] - 1] + 1] - fLoc[snpSection[snpNumber[i] - 1]]);
    if (format[0] == 4 && format[1] == 2) {
      ReadSNP42(bptr + snpLocation[snpNumber[i] - 1], numSub, dosage, p0, p1, p2);
    } else {
      infile.close();
      return 1;
    }
    for (int j = 0; j < subjectNumber.size(); ++j) {
      values(j, 4 * i) = dosage[subjectNumber[j] - 1];
      values(j, 4 * i + 1) = p0[subjectNumber[j] - 1];
      values(j, 4 * i + 2) = p1[subjectNumber[j] - 1];
      values(j, 4 * i + 3) = p2[subjectNumber[j] - 1];
    }
  }
  infile.close();
  return 0;
}