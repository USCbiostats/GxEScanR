#include <vector>
#include <iomanip>
#include "RcppArmadillo.h"

// [[Rcpp::export]]
int OpenGEOutFile(std::string &filename) {
  std::ofstream outfile;
  
  outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);
  if (!outfile.good())
    return 1;
  outfile << "SNP\tChromosome\tLocation\tReference\tAlternate\tSubjects\tCases\t";
  outfile << "betaGELin\tchiSqGELin\tbetaCaseLin\tChiSqCaseLin\tbetaControlLin\tchiSqControlLin";
  outfile << "\tbetaGELog\tchiSqGELog\tbetaCaseLog\tChiSqCaseLog\tbetaControlLog\tchiSqControlLog" << std::endl;
  outfile.close();
  return 0;
}

// [[Rcpp::export]]
int AppendGEResults(std::string &filename, Rcpp::StringVector &snpID, Rcpp::StringVector &chromosome,
                     Rcpp::IntegerVector &location, Rcpp::StringVector &refAllele, Rcpp::StringVector &altAllele,
                     int numSub, int numCases, arma::mat &logLike, arma::mat &estimates, int length,
                     double sigmaE, double sigmaEcase, double sigmaEcontrol) {
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
    outfile << estimates(i, 0) / sigmaE << '\t' << logLike(i, 0) << '\t'
            << estimates(i, 1) / sigmaEcase << '\t' << logLike(i, 1) << '\t'
            << estimates(i, 2) / sigmaEcontrol << '\t' << logLike(i, 2) << '\t'
            << estimates(i, 3) / sigmaE << '\t' << logLike(i, 3) << '\t'
            << estimates(i, 4) / sigmaEcase << '\t' << logLike(i, 4) << '\t'
            << estimates(i, 5) / sigmaEcontrol << '\t' << logLike(i, 5) << '\t'
            << std::endl;
  }
  outfile.close();
  return 0;
}
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
             arma::vec &betaT0, arma::vec &betaT, arma::vec &betaB,
             arma::mat &qL, arma::mat &qR,
             arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
             arma::vec &sigma2, arma::vec &logLike) {
  int n, p;
  
  n = xl.n_rows;

  if (solve(betaT0, rTL, qL.t() * y) == false) {
    return 1;
  }
  
  sigma2[0] = sum(square(y - xl * betaT0)) / n;
  logLike[0] = - n * (log(2. * 3.14159265) + log(sigma2[0]) + 1.) / 2.;
  if (solve(betaB, rBR, qR.t() * y) == false)
    return 1;
  if (solve(betaT, rTL, qL.t() * y - rTR * betaB) == false)
    return 1;
  sigma2[1] = sum(square(y - (xl * betaT + xr * betaB))) / n;
  logLike[1] = - n * (log(2. * 3.14159265) + log(sigma2[1]) + 1.) / 2.;
//  Rcpp::Rcout << "Linear log likelihood:\t" << logLike[0] << '\t' << logLike[1] << std::endl;
  return 0;
}

// [[Rcpp::export]]
int GELogRegStep1(arma::vec &y, arma::mat &xl, arma::mat &xr,
                  arma::vec &yp, arma::mat &xlw,
                  arma::vec &z,
                  arma::vec &beta0, arma::vec &beta0diff,
                  arma::vec &beta, arma::vec &betaT, arma::vec &betaB,
                  arma::mat &qL, arma::mat &qR,
                  arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
                  arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
                  arma::vec &w, arma::vec &wInv,
                  arma::vec &score, arma::vec &logLike) {
  abx = xl * beta0;
  expabx = exp(abx);
  expabxp1 = expabx + 1.;
  expitabx = expabx / expabxp1;
  yp = y - expitabx;
  score.subvec(0, xl.n_cols - 1) = xl.t() * yp;
  return 0;
}

// [[Rcpp::export]]
int GELogRegStep2(arma::vec &y, arma::mat &xl, arma::mat &xr,
                  arma::vec &yp, arma::mat &xlw,
                  arma::vec &z,
                  arma::vec &beta0, arma::vec &beta0diff,
                  arma::vec &beta, arma::vec &betaT, arma::vec &betaB,
                  arma::mat &qL, arma::mat &qR,
                  arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
                  arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
                  arma::vec &w, arma::vec &wInv,
                  arma::vec &score, arma::vec &logLike) {
  w = sqrt(expitabx % (1. - expitabx));
  wInv = 1. / w;
  
  return 0;
}

// [[Rcpp::export]]
int GELogRegStep3(arma::vec &y, arma::mat &xl, arma::mat &xr,
                  arma::vec &yp, arma::mat &xlw,
                  arma::vec &z,
                  arma::vec &beta0, arma::vec &beta0diff,
                  arma::vec &beta, arma::vec &betaT, arma::vec &betaB,
                  arma::mat &qL, arma::mat &qR,
                  arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
                  arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
                  arma::vec &w, arma::vec &wInv,
                  arma::vec &score, arma::vec &logLike) {
//  arma::mat res;
  
//  xlw.submat(0,0,0,0) = 1.;
//  xlw = xl.each_col() % w;
//  res = xl.each_col() % w;
//  xlw.submat(0,0,9,9) = res.submat(0,0,9,9);
//  xlw = res;
  for (int i = 0; i < xl.n_cols; ++i)
    xlw.col(i) = xl.col(i) % w;
//  xlw = 0.25 * xl;
  
  return 0;
}

// [[Rcpp::export]]
int GELogRegStep4(arma::vec &y, arma::mat &xl, arma::mat &xr,
                  arma::vec &yp, arma::mat &xlw,
                  arma::vec &z,
                  arma::vec &beta0, arma::vec &beta0diff,
                  arma::vec &beta, arma::vec &betaT, arma::vec &betaB,
                  arma::mat &qL, arma::mat &qR,
                  arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
                  arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
                  arma::vec &w, arma::vec &wInv,
                  arma::vec &score, arma::vec &logLike) {
  if (qr_econ(qL, rTL, xlw) == false)
    return 1;
  
  return 0;
}

// [[Rcpp::export]]
int GELogRegStep5(arma::vec &y, arma::mat &xl, arma::mat &xr,
                  arma::vec &yp, arma::mat &xlw,
                  arma::vec &z,
                  arma::vec &beta0, arma::vec &beta0diff,
                  arma::vec &beta, arma::vec &betaT, arma::vec &betaB,
                  arma::mat &qL, arma::mat &qR,
                  arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
                  arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
                  arma::vec &w, arma::vec &wInv,
                  arma::vec &score, arma::vec &logLike) {
  double maxBetaDiff;
  
  solve(beta0diff, rTL, qL.t() * (yp % wInv));
  maxBetaDiff = max(abs(beta0diff));
  if (maxBetaDiff > 1.) 
      beta0diff /= maxBetaDiff;
  return 0;
}

// [[Rcpp::export]]
int GELogReg(arma::vec &y, arma::mat &xl, arma::mat &xr,
             arma::vec &yp, arma::mat &xlw, arma::mat &xrw,
             arma::vec &zt, arma::vec &zb, arma::vec &k, arma::mat &t, arma::vec &h,
             arma::vec &beta0, arma::vec &beta0diff,
             arma::vec &beta, arma::vec &betaT, arma::vec &betaB,
             arma::mat &qL, arma::mat &qR,
             arma::mat &rTL, arma::mat &rTR, arma::mat &rBR,
             arma::vec &abx, arma::vec &expabx, arma::vec &expabxp1, arma::vec &expitabx,
             arma::vec &w, arma::vec &wInv,
             arma::vec &score, arma::vec &logLike) {
  int iteration;
  
  GELogRegStep1(y, xl, xr, yp, xlw, zt,
                beta0, beta0diff, beta, betaT, betaB,
                qL, qR, rTL, rTR, rBR,
                abx, expabx, expabxp1, expitabx,
                w, wInv,
                score, logLike);
//  Rcpp::Rcout << "Initial max score no E model:\t" << max(abs(score)) << std::endl;
  //if (max(abs(score)) < 1e-6)
  //  return 0;
  
  GELogRegStep2(y, xl, xr, yp, xlw, zt,
                beta0, beta0diff, beta, betaT, betaB,
                qL, qR, rTL, rTR, rBR,
                abx, expabx, expabxp1, expitabx,
                w, wInv,
                score, logLike);
  GELogRegStep3(y, xl, xr, yp, xlw, zt,
                beta0, beta0diff, beta, betaT, betaB,
                qL, qR, rTL, rTR, rBR,
                abx, expabx, expabxp1, expitabx,
                w, wInv,
                score, logLike);
  GELogRegStep4(y, xl, xr, yp, xlw, zt,
                beta0, beta0diff, beta, betaT, betaB,
                qL, qR, rTL, rTR, rBR,
                abx, expabx, expabxp1, expitabx,
                w, wInv,
                score, logLike);
  for (iteration = 0; iteration < 25 && max(abs(score)) > 1e-6; ++iteration) {
    GELogRegStep5(y, xl, xr, yp, xlw, zt,
                  beta0, beta0diff, beta, betaT, betaB,
                  qL, qR, rTL, rTR, rBR,
                  abx, expabx, expabxp1, expitabx,
                  w, wInv,
                  score, logLike);
    beta0 = beta0 + beta0diff;
    GELogRegStep1(y, xl, xr, yp, xlw, zt,
                  beta0, beta0diff, beta, betaT, betaB,
                  qL, qR, rTL, rTR, rBR,
                  abx, expabx, expabxp1, expitabx,
                  w, wInv,
                  score, logLike);
//    Rcpp::Rcout << "Max score, beta, betaDiff, yp:\t" << max(abs(score)) << '\t' << max(abs(beta0)) << '\t' << max(abs(beta0diff)) << '\t' << max(abs(yp)) << std::endl;
    if (max(abs(score)) < 1e-6)
      break;
    if (iteration < 2) {
      GELogRegStep2(y, xl, xr, yp, xlw, zt,
                    beta0, beta0diff, beta, betaT, betaB,
                    qL, qR, rTL, rTR, rBR,
                    abx, expabx, expabxp1, expitabx,
                    w, wInv,
                    score, logLike);
      GELogRegStep3(y, xl, xr, yp, xlw, zt,
                    beta0, beta0diff, beta, betaT, betaB,
                    qL, qR, rTL, rTR, rBR,
                    abx, expabx, expabxp1, expitabx,
                    w, wInv,
                    score, logLike);
      GELogRegStep4(y, xl, xr, yp, xlw, zt,
                    beta0, beta0diff, beta, betaT, betaB,
                    qL, qR, rTL, rTR, rBR,
                    abx, expabx, expabxp1, expitabx,
                    w, wInv,
                    score, logLike);
    }
  }
  if (iteration == 25)
    return 1;
  logLike[0] = sum(abx % y);
  logLike[0] -= sum(log(expabxp1));
  logLike[0] += logLike[0];
//  Rcpp::Rcout << "log likelihood:\t" << logLike[0] << std::endl;

  beta.subvec(0, beta.n_elem - 2) = beta0;
  beta(beta.n_elem - 1) = 0.;
//  Rcpp::Rcout << "Before initial max score, max score, yp:\t" << max(abs(score)) << '\t' << max(abs(yp))<< std::endl;
//  Rcpp::Rcout << "Before initial max score, max xr:\t" << max(abs(xr)) << std::endl;
//  Rcpp::Rcout << "Before initial max score, max yp:\t" << max(abs(yp)) << std::endl;
  score.subvec(score.n_elem - 1, score.n_elem - 1) = xr.t() * yp;
//  Rcpp::Rcout << "Initial max score E model:\t" << max(abs(score)) << std::endl;
  
  xrw = xr.col(0) % w;
  rTR = qL.t() * xrw;
  t = xrw - qL * rTR;
  qr_econ(qR, rBR, t);
  solve(h, rTL, rTR);

  for (iteration = 0; iteration < 25 && max(abs(score)) > 1e-6; ++iteration) {
    zt = qL.t() * (yp % wInv);
    solve(k, rTL, zt);
    zb = qR.t() * (yp % wInv);
    solve(betaB, rBR, zb);
    betaT = k - h * betaB;
    beta.subvec(0, beta.n_elem - 2) += betaT;
    beta.subvec(beta.n_elem - 1, beta.n_elem - 1) += betaB;
    
    abx = xl * beta.subvec(0, beta.n_elem - 2) + xr * beta.subvec(beta.n_elem - 1, beta.n_elem - 1);
    expabx = exp(abx);
    expabxp1 = expabx + 1.;
    expitabx = expabx / expabxp1;
    yp = y - expitabx;
    score.subvec(0, xl.n_cols - 1) = xl.t() * yp;
    score.subvec(score.n_elem - 1, score.n_elem - 1) = xr.t() * yp;
//    Rcpp::Rcout << "Max score, yp:\t" << max(abs(score)) << '\t' << max(abs(yp)) << std::endl;
  }
  if (iteration == 25)
    return 1;
  logLike[1] = sum(abx % y);
  logLike[1] -= sum(log(expabxp1));
  logLike[1] += logLike[1];
//  Rcpp::Rcout << logLike[0] << '\t' << logLike[1] << std::endl;
  
  return 0;
}

// [[Rcpp::export]]
int QREcon(arma::mat &q, arma::mat &r, arma::mat &x) {
  if (qr_econ(q, r, x) == false)
    return 1;
  return 0;
}

