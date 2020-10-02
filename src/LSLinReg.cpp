// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"

const double LOG2PIP1 = log(M_PI + M_PI) + 1.;

// [[Rcpp::export]]
Rcpp::List initlslinreg(const arma::vec &y,
                        const arma::mat &x) {
  arma::mat ql, rtl;
  arma::vec zt, k;
  double sigma2, loglike;
  
  qr_econ(ql, rtl, x);
  zt = ql.t() * y;
  solve(k, rtl, zt);
  sigma2 = arma::dot(y - x * k, y - x * k) / y.n_elem;
  loglike = -(y.n_elem * (LOG2PIP1 + log(sigma2)) / 2.);
  
  return Rcpp::List::create(Rcpp::Named("ql") = ql,
                            Rcpp::Named("rtl") = rtl,
                            Rcpp::Named("zt") = zt,
                            Rcpp::Named("k") = k,
                            Rcpp::Named("loglike") = loglike);
}

// [[Rcpp::export]]
int lslinreg(const arma::vec &y,
             const arma::mat &x,
             arma::mat &xr,
             const arma::mat &ql,
             const arma::mat &rtl,
             const arma::vec &k,
             int q,
             const Rcpp::LogicalVector &skipped,
             int skipoffset,
             int maxn,
             double minsum,
             double maxsum,
             arma::vec &loglike,
             arma::mat &beta) {
  arma::mat rtr;
  arma::mat t;
  arma::mat qr;
  arma::mat rbr;
  arma::vec zb;
  arma::vec bt;
  arma::mat h;
  arma::vec bb;
  double sigma;
  double colsum;
  double lhconst;
  int i, n;
  
  if (q < 1 || (xr.n_cols % q) != 0) {
    Rcpp::Rcerr << "Invalid number of tests to run" << std::endl;
    return 1;
  }
  
  n = xr.n_cols / q;
  if (maxn < n)
    n = maxn;
  const arma::cube xrcube(xr.memptr(), xr.n_rows, q, n, false, true);

  lhconst = -0.5 * y.n_elem*(LOG2PIP1 - log((double)y.n_elem));  
  for (i = 0; i < n; ++i) {
    if (skipped[skipoffset + i] == FALSE) {
      loglike.at(i) = NA_REAL;
      beta(i, 0) = 1.;
      continue;
    }
    colsum = sum(xrcube.slice(i).col(0));
    if (colsum < minsum || colsum > maxsum) {
      loglike.at(i) = NA_REAL;
      beta(i, 0) = 2.;
      continue;
    }
    rtr = ql.t() * xrcube.slice(i);
    t = xrcube.slice(i) - ql * rtr;
    qr_econ(qr, rbr, t);
    if (solve(h, rtl, rtr, arma::solve_opts::no_approx) == false) {
      loglike.at(i) = NA_REAL;
      beta(i, 0) = 3.;
      continue;
    }
    zb = qr.t() * y;
    if (solve(bb, rbr, zb, arma::solve_opts::no_approx) == false) {
      Rcpp::Rcout << "Solve error" << std::endl;
      beta(i, 0) = 4.;
      continue;
    }
    bt = k - h * bb;
    
    sigma = arma::norm(y - x * bt - xrcube.slice(i) * bb);
    loglike.at(i) = lhconst - y.n_elem * log(sigma);
    
    beta.row(i) = bb.t();
  }
  return 0;
}
