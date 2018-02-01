#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
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


CSubject &CSubject::operator=(const CSubject &_rhs) {
  if (this == &_rhs)
    return *this;
  m_fid = _rhs.m_fid;
  m_iid = _rhs.m_iid;
  return *this;
}

bool CSubject::operator==(const CSubject &_rhs) const {
  if (this == &_rhs)
    return true;
  return (m_fid == _rhs.m_fid && m_iid == _rhs.m_iid);
}
