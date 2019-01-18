#include <vector>
#include <iomanip>
#include "RcppArmadillo.h"

// [[Rcpp::export]]
int OpenGxEOutFile(std::string &filename) {
  std::ofstream outfile;
  
  outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
  if (!outfile.good())
    return 1;
  outfile << "SNP\tChromosome\tLocation\tReference\tAlternate\tSubjects\tCases\tbetaG\tchiSqG\tbetaGxE\tchiSqGxE\tchi2df" << std::endl;
  outfile.close();
  return 0;
}

// [[Rcpp::export]]
int AppendGxEResults(std::string &filename, Rcpp::StringVector &snpID, Rcpp::StringVector &chromosome,
                     Rcpp::IntegerVector &location, Rcpp::StringVector &refAllele, Rcpp::StringVector &altAllele,
                     int numSub, int numCases, arma::mat &logLike, arma::mat &estimates, int length, double sigmaE) {
  std::ofstream outfile;
  int i;
  
  outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::app);
  if (!outfile.good())
    return 1;

  for (i = 0; i < length; ++i) {
    if (logLike(i, 0) != logLike(i,0))
      continue;
    outfile << snpID[i] << '\t'
            << chromosome[i] << '\t'
            << location[i] << '\t'
            << refAllele[i] << '\t'
            << altAllele[i] << '\t';
    outfile << numSub << '\t' << numCases << '\t';
    outfile << estimates(i, 0) << '\t' << logLike(i, 0) << '\t'
            << estimates(i, 1) / sigmaE << '\t' << logLike(i, 1) << '\t'
            << logLike(i, 2) << std::endl;
  }
  outfile.close();
  return 0;
}
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
  if (qr_econ(ql, rtl, xlw) == false) {
    Rcpp::Rcerr << "Initialization QR decomposition failure" << std::endl;
    return 1;
  }
  yp = y - expitabx;
  zt = ql.t() * (yp % wInv);
  if (solve(k, rtl, zt) == false) {
    Rcpp::Rcerr << "Initialization solve failure for k" << std::endl;
    Rcpp::Rcerr << "k length\t" << k.n_elem << std::endl;
    Rcpp::Rcerr << "rtl\n" << rtl << std::endl;
    Rcpp::Rcerr << "zt\n" << zt << std::endl;
    return 1;
  }
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
  if (qr_econ(qr, rbr, t) == false) {
    Rcpp::Rcerr << "Maximization QR failure" << std::endl;
    Rcpp::Rcerr << "Q dimentions\t" << qr.n_rows << "\tx\t" << qr.n_cols << std::endl;
    Rcpp::Rcerr << "R dimentions\t" << rbr.n_rows << "\tx\t" << rbr.n_cols << std::endl;
    Rcpp::Rcerr << "T\n" << t << std::endl;
    return 1;
  }
  if (solve(h, rtl, rtr) == false) {
    Rcpp::Rcerr << "Maximization solve failure for h" << std::endl;
    Rcpp::Rcerr << "h dimentions\t" << h.n_rows << "\tx\t" << h.n_cols << std::endl;
    Rcpp::Rcerr << "rtl\n" << rtl << std::endl;
    Rcpp::Rcerr << "rtr\n" << rtr << std::endl;
    return 1;
  }
  
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
    if (solve(bb, rbr, zb) == false) {
      Rcpp::Rcerr << "Maximization solve failure for bb" << std::endl;
      Rcpp::Rcerr << "bb length\t" << bb.n_elem << std::endl;
      Rcpp::Rcerr << "rbr\n" << rbr << std::endl;
      Rcpp::Rcerr << "zb\n" << zb << std::endl;
      return 1;
    }
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
    if (solve(k, rtl, zt) == false) {
      Rcpp::Rcerr << "Maximization solve failure for k" << std::endl;
      Rcpp::Rcerr << "k length\t" << k.n_elem << std::endl;
      Rcpp::Rcerr << "rtl\n" << rtl << std::endl;
      Rcpp::Rcerr << "zt\n" << zt << std::endl;
      return 1;
    }
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
int ScanGenes(int n, int p, arma::vec &y, arma::mat &xl, arma::mat &xr, int numSNPs, double minMAF,
              arma::vec &beta0, arma::vec &score0, arma::vec &w, arma::vec &wInv, arma::vec &yp0,
              arma::vec &zt0, arma::vec &k0, arma::mat &ql, arma::mat &rtl,
              arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
              arma::vec &yp, arma::vec &zt, arma::vec &k, arma::vec &bt,
              arma::mat &xrw1, arma::vec &beta1, arma::vec &score1, arma::vec &zb1, arma::vec &bb1,
              arma::mat &h1, arma::mat &rtr1, arma::mat &t1, arma::mat &qr1, arma::mat &rbr1, arma::vec &logLikelihood1,
              arma::mat &xrw2, arma::vec &beta2, arma::vec &score2, arma::vec &zb2, arma::vec &bb2,
              arma::mat &h2, arma::mat &rtr2, arma::mat &t2, arma::mat &qr2, arma::mat &rbr2, arma::vec &logLikelihood2,
              arma::mat &xr1, arma::mat &xr2, arma::vec logLikelihood0, arma::mat &logLikelihoods, arma::mat &estimates) {
  int i;
  double maf;
  
//  Rcpp::Rcout << "Entering" << std::endl;

  for (i = 0; i < numSNPs; ++i) {
//    Rcpp::Rcout << i << '\t';
    xr1.submat(0, 0, n - 1, 0) = xr.submat(0, 4*i, n - 1, 4*i);
    maf = mean(xr1.col(0)) / 2.;
    if (maf < minMAF || (1. - maf) < minMAF) {
//      Rcpp::Rcout << "Skipping\t" << i << "\tMAF\t" << maf << '\t' << minMAF << std::endl;
      logLikelihoods(i,0) = NA_REAL;
      logLikelihoods(i,1) = NA_REAL;
      logLikelihoods(i,2) = NA_REAL;
      estimates(i, 0) = NA_REAL;
      estimates(i, 1) = NA_REAL;
      continue;
    }
    
    xr2.submat(0, 0, n - 1, 0) = xr.submat(0, 4*i, n - 1, 4*i);
    xr2.submat(0, 1, n - 1, 1) = xr.submat(0, 4*i, n - 1, 4*i) % xl.submat(0, xl.n_cols - 1, n - 1, xl.n_cols - 1);
    
    if (FitLRMod(n, p, y, xl, xr1,
                 beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl,
                 abx, expabx, expabxp1, expitabx, yp, zt, k, bt,
                 xrw1, beta1, score1, zb1, bb1, h1, rtr1, t1, qr1, rbr1, logLikelihood1) != 0)
      return i + 1;
    if (FitLRMod(n, p, y, xl, xr2,
                 beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl,
                 abx, expabx, expabxp1, expitabx, yp, zt, k, bt,
                 xrw2, beta2, score2, zb2, bb2, h2, rtr2, t2, qr2, rbr2, logLikelihood2) != 0)
      return i + 1;
    logLikelihoods(i, 0) = 2 * (logLikelihood1[0] - logLikelihood0[0]);
    estimates(i, 0) = beta1(p);
    logLikelihoods(i, 1) = 2 * (logLikelihood2[0] - logLikelihood1[0]);
    estimates(i, 1) = beta2(p + 1);
    logLikelihoods(i, 2) = 2 * (logLikelihood2[0] - logLikelihood0[0]);
  }  
//  Rcpp::Rcout << "Exiting" << std::endl;
  return 0;
}
