#include <RcppArmadillo.h>

int WriteOutput(const std::string &filename, const std::string &outstring);

// [[Rcpp::export]]
int InitializeLargeScaleLinReg(arma::vec &y, arma::mat &xl,
                               arma::vec &bt,
                               arma::mat &ql, arma::mat &rtl,
                               arma::vec &zt, arma::vec &k,
                               arma::vec &logLikelihood) {
  double sigma2;

  if (qr_econ(ql, rtl, xl) == false) {
    Rcpp::Rcerr << "Initialization QR decomposition failure" << std::endl;
    return 1;
  }
  zt = ql.t() * y;
  if (solve(k, rtl, zt) == false) {
    Rcpp::Rcerr << "Initialization solve failure for k" << std::endl;
    return 1;
  }
  sigma2 = arma::dot(y - xl * bt, y - xl * bt);
  sigma2 /= y.n_elem;
  logLikelihood[0] = -(log(M_PI + M_PI) + log(sigma2) + 1) / 2.;
  logLikelihood[0] *= y.n_elem;
  return 0;
}

// [[Rcpp::export]]
int ScanContinuousE(int n, int p,
                    arma::vec &y, arma::mat &xl, arma::mat &xr,
                    Rcpp::StringVector &snpID,
                    int numSNPs, double minMAF,
                    arma::mat &ql, arma::mat &rtl,
                    arma::vec &k, arma::vec &bt,
                    arma::vec &zb, arma::vec &bb,
                    arma::mat &h, arma::mat &rtr, arma::mat &t,
                    arma::mat &qr, arma::mat &rbr,
                    arma::vec &logLikelihood,
                    arma::mat &xr1,
                    arma::mat &lrTests,
                    arma::mat &estimates,
                    int testID,
                    Rcpp::StringVector &skipOut, Rcpp::StringVector &modelName) {
  int i = 0;
  double maf;
  double sigma2;
  std::string currentSNPid;
  std::ostringstream outstring;
  std::string skipFile;
  std::string modelNameC;
  
  skipFile = skipOut[0];
  modelNameC = modelName[0];
  try {
    for (i = 0; i < numSNPs; ++i) {
      lrTests(i, testID) = NA_REAL;
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
      rtr = ql.t() * xr1;
      t = xr1 - ql * rtr;
      if (qr_econ(qr, rbr, t) == false) {
        Rcpp::Rcerr << "Error in QR decomposition in " << modelNameC << " for " << snpID[i] << std::endl;
        continue;
      }
      if (solve(h, rtl, rtr) == false) {
        Rcpp::Rcerr << "Error in solve in " << modelNameC << " for " << snpID[i] << std::endl;
        continue;
      }
      zb = qr.t() * y;
      if (solve(bb, rbr, zb) == false) {
        Rcpp::Rcerr << "Error in solve in " << modelNameC << " for " << snpID[i] << std::endl;
        continue;
      }
      bt = k - h * bb;
      
      sigma2 = arma::dot(y - xl * bt - xr1 * bb, y - xl * bt - xr1 * bb);
      sigma2 /= n;
      logLikelihood[1] = -(log(M_PI + M_PI) + log(sigma2) + 1) / 2.;
      logLikelihood[1] *= n;
      lrTests(i, testID) = 2 * (logLikelihood[1] - logLikelihood[0]);
      estimates(i, testID) = bb[0];
    } 
  } catch (...) {
    Rcpp::Rcerr << snpID[i] << '\t' << "Error caught in ScanContinuousE" << std::endl;
    return 1;
  }
  
  return WriteOutput(skipFile, outstring.str());
}
