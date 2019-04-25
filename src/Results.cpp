#include <vector>
#include <iomanip>
#include "RcppArmadillo.h"

int WriteOutput(const std::string &filename, const std::string &outstring) {
  std::ofstream outfile;
  
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
  std::string outstring = "SNP\tChromosome\tLocation\tReference\tAlternate\tSubjects\tCases\tbetaG\tchiSqG\tbetaGxE\tchiSqGxE\tchi2df\n";
  
  if (filename != "stdout" && filename != "stderr") {
    outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
    if (!outfile.good()) {
      Rcpp::Rcerr << "Error opening output file" << std::endl;
      return 1;
    }
    outfile.close();
  }

  return WriteOutput(filename, outstring);
}

// [[Rcpp::export]]
int AppendGxEResults(std::string &filename, Rcpp::StringVector &snpID, Rcpp::StringVector &chromosome,
                     Rcpp::IntegerVector &location, Rcpp::StringVector &refAllele, Rcpp::StringVector &altAllele,
                     int numSub, int numCases, arma::mat &logLike, arma::mat &estimates, int length, double sigmaE) {
  std::ofstream outfile;
  std::ostringstream outstring;
  int i;
  
  for (i = 0; i < length; ++i) {
    if (logLike(i, 0) != logLike(i,0))
      continue;
    outstring << snpID[i] << '\t'
              << chromosome[i] << '\t'
              << location[i] << '\t'
              << refAllele[i] << '\t'
              << altAllele[i] << '\t';
    outstring << numSub << '\t' << numCases << '\t';
    outstring << estimates(i, 0) << '\t' << logLike(i, 0) << '\t';
    if (logLike(i, 1) != logLike(i, 1)) {
      outstring << "NA\tNA\tNA\tNA";
    } else {
      outstring << estimates(i, 1) / sigmaE << '\t' << logLike(i, 1) << '\t'
                << logLike(i, 2);
    }
    outstring << std::endl;
  }
  
  return WriteOutput(filename, outstring.str());
}