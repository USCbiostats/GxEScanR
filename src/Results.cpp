#include <vector>
#include <iomanip>
#include "RcppArmadillo.h"

int WriteOutput(const std::string &filename, const std::string &outstring) {
  std::ofstream outfile;
  
  if (filename == "")
    return 0;
  
  if (filename == "stdout") {
    Rcpp::Rcout << outstring;
    return 0;
  }
  if (filename == "stderr") {
    Rcpp::Rcerr << outstring;
    return 0;
  }
  
  outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::app);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Error opening output file" << std::endl;
    return 1;
  }
  outfile << outstring;
  outfile.close();
  
  return 0;
}
// [[Rcpp::export]]
int OpenGxEOutFile(std::string &filename) {
  std::ofstream outfile;
  std::string outstring1 = "SNP\tChromosome\tLocation\tReference\tAlternate\tSubjects\tCases\tbetaG\tchiSqG\tbetaGxE\tchiSqGxE\tchi2df\t";
  std::string outstring2 = "betaGE\tchiSqGE\tbetaCase\tchiSqCase\tbetaControl\tchiSqControl\n";
  
  if (filename != "stdout" && filename != "stderr") {
    outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
    if (!outfile.good()) {
      Rcpp::Rcerr << "Error opening output file" << std::endl;
      return 1;
    }
    outfile.close();
  }

  if (WriteOutput(filename, outstring1))
    return 1;
  return WriteOutput(filename, outstring2);
}

// [[Rcpp::export]]
int AppendGxEScanResults(std::string &filename, Rcpp::StringVector &snpID, Rcpp::StringVector &chromosome,
                         Rcpp::IntegerVector &location, Rcpp::StringVector &refAllele, Rcpp::StringVector &altAllele,
                         int numSub, int numCases, arma::mat &lrTest, arma::mat &estimates, int length, double sigmaE) {
  std::ofstream outfile;
  std::ostringstream outstring;
  int i;
  
  for (i = 0; i < length; ++i) {
    if (lrTest(i, 0) != lrTest(i,0) && lrTest(i, 3) != lrTest(i, 3))
      continue;
    outstring << snpID[i] << '\t'
              << chromosome[i] << '\t'
              << location[i] << '\t'
              << refAllele[i] << '\t'
              << altAllele[i] << '\t';
    outstring << numSub << '\t' << numCases << '\t';
    if (lrTest(i, 0) != lrTest(i, 0)) {
      outstring << "NA\tNA\tNA\tNA\tNA\tNA\t";
    } else {
      outstring << estimates(i, 0) << '\t' << lrTest(i, 0) << '\t';
      if (lrTest(i, 1) != lrTest(i, 1)) {
        outstring << "NA\tNA\tNA\t";
      } else {
        outstring << estimates(i, 1) << '\t' << lrTest(i, 1) << '\t'
                  << lrTest(i, 2) << '\t';
      }
    }
    if (lrTest(i, 3) != lrTest(i, 3)) {
      outstring << "NA\tNA\tNA\tNA\tNA\tNA";
    } else {
      outstring << estimates(i, 3) << '\t' << lrTest(i, 3) << '\t';
      if (lrTest(i, 4) != lrTest(i, 4)) {
        outstring << "NA\tNA\tNA\tNA";
      } else {
        outstring << estimates(i, 4) << '\t' << lrTest(i, 4) << '\t';
        if (lrTest(i, 5) != lrTest(i,5))
          outstring << "NA\tNA";
        else
          outstring << estimates(i, 5) << '\t' << lrTest(i, 5);
      }
    }
    outstring << std::endl;
  }
  
  return WriteOutput(filename, outstring.str());
}