#include <Rcpp.h>
#include "Subject.h"

int TestSubjectData(Rcpp::DataFrame &subjectData) {
  Rcpp::Environment env = Rcpp::Environment::namespace_env("GxEScanR");
  Rcpp::Function testSubjectFunction = env["TestSubjectDataR"];
  Rcpp::List result = testSubjectFunction(subjectData);

  Rcpp::LogicalVector success = Rcpp::as<Rcpp::LogicalVector>(result["success"]);
  if (!success[0])
    return 1;
  return 0;
}
