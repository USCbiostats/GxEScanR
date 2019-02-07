#include <vector>
#include <iomanip>
#include "RcppArmadillo.h"

// [[Rcpp::export]]
int IntializeGELinReg(arma::mat &xl, arma::mat &xr, arma::mat &ql, arma::mat &qr,
                      arma::mat &rtl, arma::mat &rtr, arma::mat &rbr) {
  arma::mat t;
  if (qr_econ(ql, rtl, xl) == false)
    return 1;
  rtr = ql.t() * xr;
  t = xr - ql * rtr;
  if (qr_econ(qr, rbr, t) == false)
    return 1;
  return 0;
}

// [[Rcpp::export]]
int GELinReg(arma::vec &y, arma::mat &xl, arma::mat &xr,
             arma::vec &betaT, arma::vec &betaB,
             arma::mat &qL, arma::mat &qR,
             arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
             arma::vec &sigma2, arma::vec &logLike) {
  int n, p;
  
  n = xl.n_rows;

  if (solve(betaT, rTL, qL.t() * y) == false) {
    return 1;
  }
  
  sigma2[0] = sum(square(y - xl * betaT)) / n;
  logLike[0] = - n * (log(2. * 3.14159265) + log(sigma2[0]) + 1.) / 2.;
  if (solve(betaB, rBR, qR.t() * y) == false)
    return 1;
  if (solve(betaT, rTL, qL.t() * y - rTR * betaB) == false)
    return 1;
  sigma2[1] = sum(square(y - (xl * betaT + xr * betaB))) / n;
  logLike[1] = - n * (log(2. * 3.14159265) + log(sigma2[1]) + 1.) / 2.;
  return 0;
}

// [[Rcpp::export]]
int GELogReg(arma::vec &y, arma::mat &xl, arma::mat &xr,
             arma::vec &yp, arma::mat &xlw,
             arma::vec &z,
             arma::vec &beta, arma::vec &betaT, arma::vec &betaB,
             arma::mat &qL, arma::mat &qR,
             arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
             arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
             arma::vec &w, arma::vec &wInv,
             arma::vec &score, arma::vec &logLike) {
  int n, p, q;
  int iterations;
  
  n = xl.n_rows;
  p = xl.n_cols;
  q = xr.n_cols;
  
  abx = xl * beta.subvec(0, p - 1);
  expabx = exp(abx);
  expabxp1 = expabx + 1.;
  expitabx = expabx / expabxp1;
  
  w = sqrt(expitabx % (1. - expabx));
  wInv = 1. / w;
  
//  xlw.submat(0, 0, 0, 0) = 5.;
//  if (p > 0)
//    return 0;
  
  Rcpp::Rcout << w.n_elem << std::endl;
  Rcpp::Rcout << xl.n_rows << '\t' << xl.n_cols << std::endl;
  Rcpp::Rcout << xlw.n_rows << '\t' << xlw.n_cols << std::endl;
  xlw.submat(0, 0, n - 1, p - 1) = xl.each_col() % w;
  Rcpp::Rcout << xlw.n_rows << '\t' << xlw.n_cols << std::endl;
  Rcpp::Rcout << xlw.submat(0, 0, 9, 9) << std::endl;
//  if (p > 0)
//    return 0;
  if (qr_econ(qL, rTL, xlw.submat(0, 0, n - 1, p - 1)) == false)
    return 1;
  if (p > 0)
    return 0;
  
  iterations = 0;
  while (iterations < 25) {
    yp = y - expitabx;
    z = qL.t() * (yp % w);
    
    score = xl.t() * yp;
    Rcpp::Rcout << "Max Score:\t" << max(abs(score)) << std::endl;
    if (max(abs(score)) < 1e-6)
      break;
    
    ++iterations;
    solve(betaT, rTL, z);
    beta.subvec(0, xl.n_cols - 1) += betaT;

    abx = xl * betaT;
    expabx = exp(abx);
    expabxp1 = expabx + 1.;
    expitabx = expabx / expabxp1;
  }
  if (iterations == 25)
    return 1;
  
  return 0;
}
