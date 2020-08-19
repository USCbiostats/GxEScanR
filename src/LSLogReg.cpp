#include "RcppArmadillo.h"

// [[Rcpp::export]]
Rcpp::List initlslogreg(const arma::vec &y,
                        const arma::mat &xl,
                        const arma::vec &beta) {
  arma::vec abx;
  arma::vec expabx;
  arma::vec expabxp1;
  arma::vec expitabx;
  arma::vec w;
  arma::vec winv;
  arma::mat xlw;
  arma::mat ql, rtl;
  arma::vec yp;
  arma::vec z;
  arma::vec k;
  arma::vec score;
  double loglike;
  
  abx = xl * beta;  
  expabx = exp(abx);
  expabxp1 = expabx + 1.;
  expitabx = expabx / expabxp1;
  w = sqrt(expitabx % (1. - expitabx));
  winv = 1. / w;
  xlw = xl.each_col() % w;
  qr_econ(ql, rtl, xlw);
  yp = y - expitabx;
  z = ql.t() * (yp % winv);
  if (solve(k, rtl, z, arma::solve_opts::no_approx) == false)  {
    Rcpp::Rcerr << "Solve error" << std::endl;
    Rcpp::stop("Error initializing large scale logistic regression");
  }
  score = xl.t() * yp;
  loglike = sum(abx % y);
  loglike -= sum(log(expabxp1));
  
  return Rcpp::List::create(Rcpp::Named("yp") = yp,
                            Rcpp::Named("ql") = ql,
                            Rcpp::Named("rtl") = rtl,
                            Rcpp::Named("k") = k,
                            Rcpp::Named("w") = w,
                            Rcpp::Named("winv") = winv,
                            Rcpp::Named("loglike") = loglike);
}


// [[Rcpp::export]]
int lslogreg(arma::vec &y,
             arma::mat &xl,
             arma::mat &xr,
             arma::vec &beta0,
             arma::vec &yp0,
             arma::mat &ql,
             arma::mat &rtl,
             arma::mat &k0,
             arma::vec &w,
             arma::vec &winv,
             int q,
             const Rcpp::LogicalVector &skipped,
             int skipoffset,
             int maxn,
             double minsum,
             double maxsum,
             arma::vec &loglike,
             arma::mat &beta) {
  
  const int maxiterations = 30;
  int iterations;
  double colsum;
  int i;
  int n;
  int ncov;
  arma::mat xrw;
  arma::mat rtr;
  arma::mat t;
  arma::mat qr;
  arma::mat rbr;
  arma::mat h;
  arma::vec betat;
  arma::vec betab;
  arma::vec yp;
  arma::vec scoret;
  arma::vec scoreb;
  arma::vec k;
  arma::vec zt;
  arma::vec zb;
  arma::vec bt;
  arma::vec bb;
  arma::vec abx;
  arma::vec expabx;
  arma::vec expabxp1;
  arma::vec expitabx;
  
  if (q < 1 || (xr.n_cols % q) != 0) {
    Rcpp::Rcerr << "Invalid number of tests to run" << std::endl;
    return 1;
  }

  n = xr.n_cols / q;
  if (maxn < n)
    n = maxn;
  ncov = xl.n_cols + q;

  const arma::cube xrcube(xr.memptr(), xr.n_rows, q, n, false, true);
  beta.zeros(n, q);
  loglike.zeros(n);

  for (i = 0; i < n; ++i) {
    if (skipped[skipoffset + i] == FALSE) {
      loglike(i) = NA_REAL;
      beta(i, 0) = 1.;
      continue;
    }
    colsum = sum(xrcube.slice(i).col(0));
    if (colsum < minsum || colsum > maxsum) {
      loglike(i) = NA_REAL;
      beta(i, 0) = 2.;
      continue;
    }
    
    xrw = xrcube.slice(i).each_col() % w;
    rtr = ql.t() * xrw;
    t = xrw - ql * rtr;
    qr_econ(qr, rbr, t);
    if (solve(h, rtl, rtr, arma::solve_opts::no_approx) == false) {
      loglike.at(i) = NA_REAL;
      beta(i, 0) = 3.;
      continue;
    }

    betat = beta0;
    betab.zeros(q);

    scoret.zeros(ncov);
    scoreb = xrcube.slice(i).t() * yp0;
    if (max(abs(scoreb)) < 0.000001) {
      abx = xl * betat + xrcube.slice(i) * betab;  
      expabx = exp(abx);
      expabxp1 = expabx + 1.;
      expitabx = expabx / expabxp1;
      loglike(i) = sum(abx % y) - sum(log(expabxp1));
      beta.row(i) = betab.t();
      continue;
    }
    beta.row(i) = scoreb.t();
    
    betat = beta0;
    betab.zeros(q);
    yp = yp0 % winv;
    k = k0;

    iterations = 0;
    while (iterations < maxiterations) {
      zb = qr.t() * yp;
      if (solve(bb, rbr, zb, arma::solve_opts::no_approx) == false) {
        loglike(i) = NA_REAL;
        beta(i, 0) = 4.;
        continue;
      }
      bt = k - h * bb;
      betat += bt;
      betab += bb;
      
      abx = xl * betat + xrcube.slice(i) * betab;  
      expabx = exp(abx);
      expabxp1 = expabx + 1.;
      expitabx = expabx / expabxp1;
      yp = y - expitabx;
      scoret = xl.t() * yp;
      scoreb = xrcube.slice(i).t() * yp;
      
      if (max(abs(scoret)) < 1e-6 && max(abs(scoreb)) < 1e-6) {
        loglike(i) = sum(abx % y) - sum(log(expabxp1));
        beta.row(i) = betab.t();
        break;
      }
      
      yp %= winv;
      zt = ql.t() * yp;
      if (solve(k, rtl, zt, arma::solve_opts::no_approx) == false) {
        loglike(i) = NA_REAL;
        beta(i, 0) = 5.;
        continue;
      }
      ++iterations;
    }
    if (iterations == maxiterations) {
      loglike(i) = NA_REAL;
      beta(i, 0) = 6.;
    }
  }
  return 0;
}
