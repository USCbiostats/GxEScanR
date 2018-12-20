// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <iostream>
#include <fstream>
#include <RcppArmadillo.h>

int FitCovariateOnlyModel(const arma::vec &y, const arma::mat &x, arma::vec &beta) {
  arma::vec sc;
  arma::vec abx, expabx, expabxp1, expitabx;
  arma::vec score;
  arma::vec w;
  arma::mat wx;
  arma::vec yp, ypw;
  arma::mat q, r;
  arma::vec betadiff;
  double scoreNorm0, scoreNorm1;

  sc = x.t() * y;
  abx = x * beta;
  expabx = exp(abx);
  expabxp1 = expabx + 1.;
  expitabx = expabx / expabxp1;
  score = sc - x.t() * expitabx;
  
  if (max(abs(score)) < 1e-6)
    return 0;
//  scoreNorm0 = norm(score);
  scoreNorm0 = max(abs(score));
  Rcpp::Rcout << scoreNorm0 << std::endl;

  w = sqrt(expitabx % (1. - expitabx));
  wx = x.each_col() % w;
  yp = (y - expitabx) / w;
  qr_econ(q, r, wx);

  for (int i = 1; i < 11; ++i) {
    ypw = q.t() * yp;
    betadiff = arma::solve(r, ypw);
    beta += betadiff;
    if (i == 1)
      Rcpp::Rcout << "Beta diff\n" << betadiff.subvec(0,4) << std::endl;

    abx = x * beta;
    expabx = exp(abx);
    expabxp1 = expabx + 1.;
    expitabx = expabx / expabxp1;
    score = sc - x.t() * expitabx;
    if (max(abs(score)) < 1e-6)
      return 0;
    //  scoreNorm1 = norm(score);
    scoreNorm1 = max(abs(score));
    scoreNorm0 = scoreNorm1;
    Rcpp::Rcout << scoreNorm0 << std::endl;
    
    yp = (y - expitabx) / w;
  }
//  Rcpp::Rcout << std::endl << beta << std::endl;
  return 0;
}

int FitCovariateOnlyModel2(const arma::vec &y, const arma::mat &x, arma::vec &beta) {
  arma::vec sc;
  arma::vec abx, expabx, expabxp1, expitabx;
  arma::mat xt;
  arma::vec score;
  arma::mat info;
  arma::mat q, r;
  arma::vec betadiff, qbeta;
  double scoreNorm0, scoreNorm1;
  arma::vec betaNew;
  arma::vec expitNew;
  arma::vec p;
  arma::mat wx;
  
  sc = x.t() * y;
  abx = x * beta;
  expabx = exp(abx);
  expabxp1 = expabx + 1.;
  expitabx = expabx / expabxp1;
  
  score = sc - x.t() * expitabx;
  if (max(abs(score)) < 1e-6)
    return 0;
  //  scoreNorm0 = norm(score);
  scoreNorm0 = max(abs(score));
  Rcpp::Rcout << scoreNorm0 << std::endl;

  xt = x.t();
  p = expitabx % (1. - expitabx);
  wx = x.each_col() % p;
  info = xt * wx;
//  Rcpp::Rcout << info.submat(0, 0, 4, 4) << std::endl;

  qr_econ(q, r, info);
//  Rcpp::Rcout << "Q\n" << q.submat(0,0,4, 4); 
//  Rcpp::Rcout << "R\n" << r.submat(0,0,4, 4); 
  
  for (int i = 1; i < 20; ++i) {
    qbeta = q.t() * score;
//    Rcpp::Rcout << "QtBeta\n" << qbeta.subvec(0,4); 
    betadiff = arma::solve(r, qbeta);
    if (i == 1)
      Rcpp::Rcout << "Beta diff\n" << betadiff.subvec(0,4) << std::endl;
    betaNew = beta + betadiff;
    //    Rcpp::Rcout << "After beta" << std::endl;
    
    abx = x * betaNew;
    expabx = exp(abx);
    expabxp1 = expabx + 1.;
    expitNew = expabx / expabxp1;
    
    score = sc - x.t() * expitNew;
    //  scoreNorm1 = norm(score);
    scoreNorm1 = max(abs(score));
    Rcpp::Rcout << scoreNorm1 << std::endl;
    if (max(abs(score)) < 1e-6)
      break;
/*    
    if (scoreNorm1 > scoreNorm0) {
      Rcpp::Rcout << "Updating" << std::endl;
      wx = x.each_col() % w;
      qr_econ(q, r, wx);
      
      ypw = q.t() * yp;
      betadiff = arma::solve(r, ypw);
      betaNew = beta + betadiff;
      
      abx = x * betaNew;
      expabx = exp(abx);
      expabxp1 = expabx + 1.;
      expitNew = expabx / expabxp1;
      
      score = sc - x.t() * expitNew;
      //  scoreNorm1 = norm(score);
      scoreNorm1 = max(abs(score));
      Rcpp::Rcout << "Updated\t" << scoreNorm1 << std::endl;
      if (max(abs(score)) < 1e-6)
        break;
    }
*/
    scoreNorm0 = scoreNorm1;
    
    beta = betaNew;
    expitabx = expitNew;
  }
  Rcpp::Rcout << std::endl << beta.subvec(0,4) << std::endl;

  return 0;
}

// [[Rcpp::export]]
int GxEScanC(arma::vec &phenotype, arma::mat &covariates, arma::vec &mu, arma::vec &sigma2, arma::vec &beta0, arma::icolvec &index,
             Rcpp::List &geneticData, std::string &outFilename, std::string &skippedFilename, double minmaf) {
  arma::vec beta;
  
  if (phenotype.size() == 0) {
    Rcpp::Rcerr << "No subjects have complete data" << std::endl;
    return 1;
  }
  beta = beta0;
  beta[beta.n_elem - 1] = 0.;
  beta.zeros(covariates.n_cols);
  FitCovariateOnlyModel(phenotype, covariates, beta);
  beta.zeros(covariates.n_cols);
  FitCovariateOnlyModel2(phenotype, covariates, beta);
  return 0;
}