#include <Rcpp.h>
#include "Subject.h"

int TestSubjectData() {
  Rcpp::Environment env("package:GxEScanR");
  Rcpp::Rcout << "Environment" << std::endl;
  Rcpp::Function testSubjectFunction = env["TestSubjectDataR"];
  Rcpp::Rcout << "Function" << std::endl;
  Rcpp::IntegerVector result = testSubjectFunction();
  
  if (result[0] == 0)
    Rcpp::Rcout << "Subject data good" << std::endl;
  else
    Rcpp::Rcout << "Subject data error" << std::endl;
  return result[0];
}
