#include <RcppArmadillo.h>

int WriteOutput(const std::string &filename, const std::string &outstring);

// [[Rcpp::export]]
int InitializeLargeScaleLogReg(int numRow, int numCol, arma::vec &y, arma::mat &xl,
                               arma::vec &beta, arma::vec &score, arma::vec &w, arma::vec &wInv,
                               arma::vec &yp, arma::vec &zt, arma::vec &k,
                               arma::mat &ql, arma::mat &rtl,
                               arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1,
                               arma::vec &expitabx, arma::vec &logLikelihood) {
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
    return 1;
  }
  score = xl.t() * yp;
  logLikelihood[0] = sum(abx % y);
  logLikelihood[0] -= sum(log(expabxp1));
  return 0;
}


/// [[Rcpp::export]]
int FitLRMod(int n, int p,
             arma::vec &y, arma::mat &xl, arma::mat &xr,
             std::string &snpID, std::string modelString,
             arma::vec &beta0, arma::vec &score0,
             arma::vec &w, arma::vec &wInv,
             arma::vec &yp0, arma::vec &zt0, arma::vec &k0,
             arma::mat &ql, arma::mat &rtl,
             arma::vec &abx, arma::vec &expabx,
             arma::vec &expabxp1, arma::vec &expitabx,
             arma::vec &yp, arma::vec &zt, arma::vec &k,
             arma::vec &bt, arma::mat &xrw, arma::vec &beta,
             arma::vec &score, arma::vec &zb, arma::vec &bb,
             arma::mat &h, arma::mat &rtr, arma::mat &t,
             arma::mat &qr, arma::mat &rbr, arma::vec &logLikelihood,
             std::ostream &outstream) {
  const int maxIterations = 20;
  int iterations;
  int q;
  
  logLikelihood[0] = NA_REAL;
  try {
    q = p + xr.n_cols;
    
    xrw = xr.each_col() % w;
    rtr = ql.t() * xrw;
    t = xrw - ql * rtr;
    if (qr_econ(qr, rbr, t) == false) {
      outstream << snpID << '\t' << "Maximization QR failure for " << modelString << std::endl;
//      Rcpp::Rcerr << "Q dimentions\t" << qr.n_rows << "\tx\t" << qr.n_cols << std::endl;
//      Rcpp::Rcerr << "R dimentions\t" << rbr.n_rows << "\tx\t" << rbr.n_cols << std::endl;
//      Rcpp::Rcerr << "T\n" << t << std::endl;
      return 1;
    }
    if (solve(h, rtl, rtr) == false) {
      outstream << snpID << '\t' << "Maximization solve failure for h for " << modelString << std::endl;
//      Rcpp::Rcerr << "h dimentions\t" << h.n_rows << "\tx\t" << h.n_cols << std::endl;
//      Rcpp::Rcerr << "rtl\n" << rtl << std::endl;
//      Rcpp::Rcerr << "rtr\n" << rtr << std::endl;
      return 1;
    }
    
    beta.zeros();
    beta.subvec(0, beta0.n_elem - 1) = beta0;
    
    yp = yp0;
    score.subvec(0, score0.n_elem - 1) = score0;
    score.subvec(p, q - 1) = xr.t() * yp;
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
        outstream << snpID << '\t' << "Maximization solve failure for bb for " << modelString << std::endl;
//        Rcpp::Rcerr << "bb length\t" << bb.n_elem << std::endl;
//        Rcpp::Rcerr << "rbr\n" << rbr << std::endl;
//        Rcpp::Rcerr << "zb\n" << zb << std::endl;
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
        outstream << snpID << '\t' << "Maximization solve failure for k for " << modelString << std::endl;
//        Rcpp::Rcerr << "k length\t" << k.n_elem << std::endl;
//        Rcpp::Rcerr << "rtl\n" << rtl << std::endl;
//        Rcpp::Rcerr << "zt\n" << zt << std::endl;
        return 1;
      }
      ++iterations;
    }
    if (iterations != maxIterations) {
      logLikelihood[0] = sum(abx % y);
      logLikelihood[0] -= sum(log(expabxp1));
    } else {
      outstream << snpID << '\t' << "Maximum iterations exceeded for " << modelString << std::endl;
      return 1;
    }
  } catch (...) {
    outstream << snpID << '\t' << "Error caught in FItLRMod for " << modelString << std::endl;
    return 1;
  }
  return 0;
}

// [[Rcpp::export]]
int ScanDisease(int n, int p,
                arma::vec &y, arma::mat &xl, arma::mat &xr,
                Rcpp::StringVector &snpID,
                int numSNPs, double minMAF,
                arma::vec &beta0, arma::vec &score0,
                arma::vec &w, arma::vec &wInv,
                arma::vec &yp0, arma::vec &zt0, arma::vec &k0,
                arma::mat &ql, arma::mat &rtl,
                arma::vec &abx, arma::vec &expabx,
                arma::vec &expabxp1, arma::vec &expitabx,
                arma::vec &yp, arma::vec &zt, arma::vec &k, arma::vec &bt,
                arma::mat &xrw1, arma::vec &beta1, arma::vec &score1,
                arma::vec &zb1, arma::vec &bb1,
                arma::mat &h1, arma::mat &rtr1, arma::mat &t1,
                arma::mat &qr1, arma::mat &rbr1,
                arma::vec &logLikelihood1,
                arma::mat &xrw2, arma::vec &beta2, arma::vec &score2,
                arma::vec &zb2, arma::vec &bb2,
                arma::mat &h2, arma::mat &rtr2, arma::mat &t2,
                arma::mat &qr2, arma::mat &rbr2,
                arma::vec &logLikelihood2,
                arma::mat &xr1, arma::mat &xr2,
                arma::vec logLikelihood0, arma::mat &logLikelihoods,
                arma::mat &estimates,
                Rcpp::StringVector &skipOut) {
  int i = 0;
  double maf;
  std::string currentSNPid;
  std::string modelDG = "model D|G";
  std::string modelDGxE = "model D|G,GxE";
  std::ostringstream outstring;
  std::string skipFile;
  
  skipFile = skipOut[0];
  try {
    for (i = 0; i < numSNPs; ++i) {
      logLikelihoods(i,0) = NA_REAL;
      logLikelihoods(i,1) = NA_REAL;
      logLikelihoods(i,2) = NA_REAL;
      estimates(i, 0) = NA_REAL;
      estimates(i, 1) = NA_REAL;

      xr1.submat(0, 0, n - 1, 0) = xr.submat(0, 4*i, n - 1, 4*i);
      maf = mean(xr1.col(0)) / 2.;
      if (maf < minMAF || (1. - maf) < minMAF) {
        if (skipFile != "") {
          outstring << snpID[i] << "\tskipped because maf, ";
          outstring << (maf < 0.5 ? maf : 1 - maf);
          outstring << ", is less than sampleminMaf, " << minMAF << std::endl;
        }
        continue;
      }
      
      xr2.submat(0, 0, n - 1, 0) = xr.submat(0, 4*i, n - 1, 4*i);
      xr2.submat(0, 1, n - 1, 1) = xr.submat(0, 4*i, n - 1, 4*i) % xl.submat(0, xl.n_cols - 1, n - 1, xl.n_cols - 1);
      currentSNPid = snpID[i];
      if (FitLRMod(n, p, y, xl, xr1, currentSNPid, modelDG,
                   beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl,
                   abx, expabx, expabxp1, expitabx, yp, zt, k, bt,
                   xrw1, beta1, score1, zb1, bb1, h1, rtr1, t1, qr1,
                   rbr1, logLikelihood1, outstring) == 0) {
        logLikelihoods(i, 0) = 2 * (logLikelihood1[0] - logLikelihood0[0]);
        estimates(i, 0) = beta1(p);
        if (FitLRMod(n, p, y, xl, xr2, currentSNPid, modelDGxE,
                     beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl,
                     abx, expabx, expabxp1, expitabx, yp, zt, k, bt,
                     xrw2, beta2, score2, zb2, bb2, h2, rtr2, t2, qr2,
                     rbr2, logLikelihood2, outstring) == 0) {
          logLikelihoods(i, 0) = 2 * (logLikelihood1[0] - logLikelihood0[0]);
          estimates(i, 0) = beta1(p);
          logLikelihoods(i, 1) = 2 * (logLikelihood2[0] - logLikelihood1[0]);
          estimates(i, 1) = beta2(p + 1);
          logLikelihoods(i, 2) = 2 * (logLikelihood2[0] - logLikelihood0[0]);
        }
      }
    } 
  } catch (...) {
    Rcpp::Rcerr << snpID[i] << '\t' << "Error caught in ScanGenes" << std::endl;
    return 1;
  }
  
  return WriteOutput(skipFile, outstring.str());
}

// [[Rcpp::export]]
int ScanBinaryE(int n, int p,
                arma::vec &y, arma::mat &xl, arma::mat &xr,
                Rcpp::StringVector &snpID,
                int numSNPs, double minMAF,
                arma::vec &beta0, arma::vec &score0,
                arma::vec &w, arma::vec &wInv,
                arma::vec &yp0, arma::vec &zt0, arma::vec &k0,
                arma::mat &ql, arma::mat &rtl,
                arma::vec &abx, arma::vec &expabx,
                arma::vec &expabxp1, arma::vec &expitabx,
                arma::vec &yp, arma::vec &zt, arma::vec &k, arma::vec &bt,
                arma::mat &xrw1, arma::vec &beta1, arma::vec &score1,
                arma::vec &zb1, arma::vec &bb1,
                arma::mat &h1, arma::mat &rtr1, arma::mat &t1,
                arma::mat &qr1, arma::mat &rbr1,
                arma::vec &logLikelihood1,
                arma::mat &xr1,
                arma::vec logLikelihood0, arma::mat &logLikelihoods,
                arma::mat &estimates, int testID,
                Rcpp::StringVector &skipOut, Rcpp::StringVector &modelName) {
  int i = 0;
  double maf;
  std::string currentSNPid;
  std::ostringstream outstring;
  std::string skipFile;
  std::string modelNameC;
  
  skipFile = skipOut[0];
  modelNameC = modelName[0];
  try {
    for (i = 0; i < numSNPs; ++i) {
      logLikelihoods(i, testID) = NA_REAL;
      estimates(i, testID) = NA_REAL;

      xr1.submat(0, 0, n - 1, 0) = xr.submat(0, 4*i, n - 1, 4*i);
      maf = mean(xr1.col(0)) / 2.;
      if (maf < minMAF || (1. - maf) < minMAF) {
// This is currently commented out but may be readded if
// option to only do E|G is added.
//        outstring << snpID[i] << "\tskipped because maf, ";
//        outstring << (maf < 0.5 ? maf : 1 - maf);
//        outstring << ", is less than sampleminMaf, " << minMAF << std::endl;
        continue;
      }
      
      currentSNPid = snpID[i];
      if (FitLRMod(n, p, y, xl, xr1, currentSNPid, modelNameC,
                   beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl,
                   abx, expabx, expabxp1, expitabx, yp, zt, k, bt,
                   xrw1, beta1, score1, zb1, bb1, h1, rtr1, t1, qr1,
                   rbr1, logLikelihood1, outstring) == 0) {
        logLikelihoods(i, testID) = 2 * (logLikelihood1[0] - logLikelihood0[0]);
        estimates(i, testID) = beta1(p);
      }
    } 
  } catch (...) {
    Rcpp::Rcerr << snpID[i] << '\t' << "Error caught in ScanBinaryE" << std::endl;
    return 1;
  }
  
  return WriteOutput(skipFile, outstring.str());
}
