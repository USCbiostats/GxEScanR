#include <vector>
#include <iomanip>
#include "RcppArmadillo.h"

// [[Rcpp::export]]
int InitializeLRMod(int numRow, int numCol, arma::vec &y, arma::mat &xl,
                    arma::vec &beta, arma::vec &score, arma::vec &w, arma::vec &wInv,
                    arma::vec &yp, arma::vec &zt, arma::vec &k,
                    arma::mat &ql, arma::mat &rtl,
                    arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx, arma::vec &logLikelihood) {
  arma::mat xlw;

  abx = xl * beta;  
  expabx = exp(abx);
  expabxp1 = expabx + 1.;
  expitabx = expabx / expabxp1;
  w = sqrt(expitabx % (1. - expitabx));
  wInv = 1. / w;
  xlw = xl.each_col() % w;
  qr_econ(ql, rtl, xlw);
  yp = y - expitabx;
  zt = ql.t() * (yp % wInv);
  k = solve(rtl, zt);
  score = xl.t() * yp;
  logLikelihood[0] = sum(abx % y);
  logLikelihood[0] -= sum(log(expabxp1));
  return 0;
}

// [[Rcpp::export]]
int FitLRMod(int n, int p, arma::vec &y, arma::mat &xl, arma::mat &xr,
             arma::vec &beta0, arma::vec &score0, arma::vec &w, arma::vec &wInv, arma::vec &yp0,
             arma::vec &zt0, arma::vec &k0, arma::mat &ql, arma::mat &rtl,
             arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
             arma::vec &yp, arma::vec &zt, arma::vec &k, arma::vec &bt,
             arma::mat &xrw, arma::vec &beta, arma::vec &score, arma::vec &zb, arma::vec &bb,
             arma::mat &h, arma::mat &rtr, arma::mat &t, arma::mat &qr, arma::mat &rbr, arma::vec &logLikelihood) {
  const int maxIterations = 5;
  int iterations;
  int q;
  
  q = p + xr.n_cols;
  
  xrw = xr.each_col() % w;
  rtr = ql.t() * xrw;
  t = xrw - ql * rtr;
  qr_econ(qr, rbr, t);
  h = solve(rtl, rtr);
  
  beta.zeros();
  beta.subvec(0, beta0.n_elem - 1) = beta0;

  yp = yp0;
  score.subvec(0, score0.n_elem - 1) = score0;
  score.subvec(p, q - 1) = xr.t() * yp;
//  Rcpp::Rcout << max(abs(score)) << std::endl;
  // Check for convergence
  if (max(abs(score)) < 1e-6)
    return 0;
  yp %= wInv;
  zt = zt0;
  k = k0;

  iterations = 0;
  while (iterations < maxIterations) {
    // Update beta
    zb = qr.t() * yp;
    bb = solve(rbr, zb);
    bt = k - h * bb;
    beta.subvec(0, p - 1) += bt;
    beta.subvec(p, q - 1) += bb;
    
    // Calculate the score
    abx = xl * beta.subvec(0, p - 1) + xr * beta.subvec(p, q - 1);  
    expabx = exp(abx);
    expabxp1 = expabx + 1.;
    expitabx = expabx / expabxp1;
    yp = y - expitabx;
    score.subvec(0, p - 1) = xl.t() * yp;
    score.subvec(p, q - 1) = xr.t() * yp;
//    Rcpp::Rcout << max(abs(score)) << std::endl;
    // Check for convergence
    if (max(abs(score)) < 1e-6)
      break;
    
    // Update the bits needed for the next iteration
    yp %= wInv;
    zt = ql.t() * yp;
    k = solve(rtl, zt);
    ++iterations;
  }
  logLikelihood[0] = sum(abx % y);
  logLikelihood[0] -= sum(log(expabxp1));
//  Rcpp::Rcout << std::setprecision(10) << sum(abx % y) << '\t' << std::setprecision(10) << sum(log(expabxp1)) << std::endl;
//  Rcpp::Rcout << std::setprecision(10) << logLikelihood[0] << std::endl;
//  Rcpp::Rcout << min(abs(log(expabxp1))) << '\t' << max(abs(log(expabxp1))) << std::endl;
  return 0;
}

// [[Rcpp::export]]
Rcpp::List ScanGenes(int n, int p, arma::vec &y, arma::mat &xl, arma::mat &xr, int numSNPs,
                     arma::vec &beta0, arma::vec &score0, arma::vec &w, arma::vec &wInv, arma::vec &yp0,
                     arma::vec &zt0, arma::vec &k0, arma::mat &ql, arma::mat &rtl,
                     arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
                     arma::vec &yp, arma::vec &zt, arma::vec &k, arma::vec &bt,
                     arma::mat &xrw1, arma::vec &beta1, arma::vec &score1, arma::vec &zb1, arma::vec &bb1,
                     arma::mat &h1, arma::mat &rtr1, arma::mat &t1, arma::mat &qr1, arma::mat &rbr1, arma::vec &logLikelihood1,
                     arma::mat &xrw2, arma::vec &beta2, arma::vec &score2, arma::vec &zb2, arma::vec &bb2,
                     arma::mat &h2, arma::mat &rtr2, arma::mat &t2, arma::mat &qr2, arma::mat &rbr2, arma::vec &logLikelihood2) {
  Rcpp::List retVal;
  int i;
  Rcpp::Rcout << "Entering" << std::endl;
  arma::mat xr1(n, 1, arma::fill::zeros);
  arma::mat xr2(n, 2, arma::fill::zeros);
  arma::mat logLikelihood(numSNPs, 2);

  for (i = 0; i < numSNPs; ++i) {
    xr1.submat(0, 0, n - 1, 0) = xr.submat(0, 0, n - 1, 0);
    xr2.submat(0, 0, n - 1, 0) = xr.submat(0, 0, n - 1, 0);
    xr2.submat(0, 1, n - 1, 1) = xr.submat(0, 0, n - 1, 0) % xl.submat(0, xl.n_cols - 1, n - 1, xl.n_cols - 1);
    
    FitLRMod(n, p, y, xl, xr1,
             beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl,
             abx, expabx, expabxp1, expitabx, yp, zt, k, bt,
             xrw1, beta1, score1, zb1, bb1, h1, rtr1, t1, qr1, rbr1, logLikelihood1);
    FitLRMod(n, p, y, xl, xr2,
             beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl,
             abx, expabx, expabxp1, expitabx, yp, zt, k, bt,
             xrw2, beta2, score2, zb2, bb2, h2, rtr2, t2, qr2, rbr2, logLikelihood2);
    logLikelihood(0, 0) = logLikelihood1[0];
    logLikelihood(0, 1) = logLikelihood2[0];
  }  
  
  retVal = Rcpp::List::create(Rcpp::Named("LogLikelihoods") = logLikelihood);
  return retVal;
  //  return Rcpp::List::create(Rcpp::Named("XR1") = xr1,
//                               Rcpp::Named("XR2") = xr2);
}