#include <vector>
#include <iomanip>
#include "RcppArmadillo.h"

// [[Rcpp::export]]
int OpenGxEOutFile(std::string &filename) {
  std::ofstream outfile;
  
  if (filename == "NoOutput") {
    return 0;
  }
  
  if (filename == "RTerminal") {
    Rcpp::Rcout << "SNP\tChromosome\tLocation\tReference\tAlternate\tSubjects\tCases\tbetaG\tchiSqG\tbetaGxE\tchiSqGxE\tchi2df" << std::endl;
    return 0;
  }
  
  outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
  if (!outfile.good())
    return 1;
  outfile << "SNP\tChromosome\tLocation\tReference\tAlternate\tSubjects\tCases\tbetaG\tchiSqG\tbetaGxE\tchiSqGxE\tchi2df" << std::endl;
  outfile.close();
  return 0;
}
