// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// initlslinreg
Rcpp::List initlslinreg(const arma::vec& y, const arma::mat& x);
RcppExport SEXP _GxEScanR_initlslinreg(SEXP ySEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(initlslinreg(y, x));
    return rcpp_result_gen;
END_RCPP
}
// lslinreg
int lslinreg(const arma::vec& y, const arma::mat& x, arma::mat& xr, const arma::mat& ql, const arma::mat& rtl, const arma::vec& k, int q, const Rcpp::LogicalVector& skipped, int skipoffset, int maxn, double minsum, double maxsum, arma::vec& loglike, arma::mat& beta);
RcppExport SEXP _GxEScanR_lslinreg(SEXP ySEXP, SEXP xSEXP, SEXP xrSEXP, SEXP qlSEXP, SEXP rtlSEXP, SEXP kSEXP, SEXP qSEXP, SEXP skippedSEXP, SEXP skipoffsetSEXP, SEXP maxnSEXP, SEXP minsumSEXP, SEXP maxsumSEXP, SEXP loglikeSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rtl(rtlSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type skipped(skippedSEXP);
    Rcpp::traits::input_parameter< int >::type skipoffset(skipoffsetSEXP);
    Rcpp::traits::input_parameter< int >::type maxn(maxnSEXP);
    Rcpp::traits::input_parameter< double >::type minsum(minsumSEXP);
    Rcpp::traits::input_parameter< double >::type maxsum(maxsumSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type loglike(loglikeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(lslinreg(y, x, xr, ql, rtl, k, q, skipped, skipoffset, maxn, minsum, maxsum, loglike, beta));
    return rcpp_result_gen;
END_RCPP
}
// initlslogreg
Rcpp::List initlslogreg(const arma::vec& y, const arma::mat& xl, const arma::vec& beta);
RcppExport SEXP _GxEScanR_initlslogreg(SEXP ySEXP, SEXP xlSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(initlslogreg(y, xl, beta));
    return rcpp_result_gen;
END_RCPP
}
// lslogreg
int lslogreg(arma::vec& y, arma::mat& xl, arma::mat& xr, arma::vec& beta0, arma::vec& yp0, arma::mat& ql, arma::mat& rtl, arma::mat& k0, arma::vec& w, arma::vec& winv, int q, const Rcpp::LogicalVector& skipped, int skipoffset, int maxn, double minsum, double maxsum, arma::vec& loglike, arma::mat& beta);
RcppExport SEXP _GxEScanR_lslogreg(SEXP ySEXP, SEXP xlSEXP, SEXP xrSEXP, SEXP beta0SEXP, SEXP yp0SEXP, SEXP qlSEXP, SEXP rtlSEXP, SEXP k0SEXP, SEXP wSEXP, SEXP winvSEXP, SEXP qSEXP, SEXP skippedSEXP, SEXP skipoffsetSEXP, SEXP maxnSEXP, SEXP minsumSEXP, SEXP maxsumSEXP, SEXP loglikeSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yp0(yp0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtl(rtlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type winv(winvSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type skipped(skippedSEXP);
    Rcpp::traits::input_parameter< int >::type skipoffset(skipoffsetSEXP);
    Rcpp::traits::input_parameter< int >::type maxn(maxnSEXP);
    Rcpp::traits::input_parameter< double >::type minsum(minsumSEXP);
    Rcpp::traits::input_parameter< double >::type maxsum(maxsumSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type loglike(loglikeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(lslogreg(y, xl, xr, beta0, yp0, ql, rtl, k0, w, winv, q, skipped, skipoffset, maxn, minsum, maxsum, loglike, beta));
    return rcpp_result_gen;
END_RCPP
}
// readblock
void readblock(Rcpp::StringVector& filename, const double blkloc, const double blkbytes, arma::ivec& blkbuffer);
RcppExport SEXP _GxEScanR_readblock(SEXP filenameSEXP, SEXP blklocSEXP, SEXP blkbytesSEXP, SEXP blkbufferSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const double >::type blkloc(blklocSEXP);
    Rcpp::traits::input_parameter< const double >::type blkbytes(blkbytesSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type blkbuffer(blkbufferSEXP);
    readblock(filename, blkloc, blkbytes, blkbuffer);
    return R_NilValue;
END_RCPP
}
// getdosages
void getdosages(arma::mat& dosages, arma::ivec& blkbuffer, double fileloc, arma::vec& indices, int firstsnp, int lastsnp, int base);
RcppExport SEXP _GxEScanR_getdosages(SEXP dosagesSEXP, SEXP blkbufferSEXP, SEXP filelocSEXP, SEXP indicesSEXP, SEXP firstsnpSEXP, SEXP lastsnpSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type dosages(dosagesSEXP);
    Rcpp::traits::input_parameter< arma::ivec& >::type blkbuffer(blkbufferSEXP);
    Rcpp::traits::input_parameter< double >::type fileloc(filelocSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< int >::type firstsnp(firstsnpSEXP);
    Rcpp::traits::input_parameter< int >::type lastsnp(lastsnpSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    getdosages(dosages, blkbuffer, fileloc, indices, firstsnp, lastsnp, base);
    return R_NilValue;
END_RCPP
}
// stdmat
void stdmat(arma::mat& m1, arma::mat& m2, arma::rowvec& means, arma::rowvec& stddevs);
RcppExport SEXP _GxEScanR_stdmat(SEXP m1SEXP, SEXP m2SEXP, SEXP meansSEXP, SEXP stddevsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type m2(m2SEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type means(meansSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type stddevs(stddevsSEXP);
    stdmat(m1, m2, means, stddevs);
    return R_NilValue;
END_RCPP
}
// makegxr
void makegxr(arma::mat& dest, const arma::mat& src, const arma::uvec& idx);
RcppExport SEXP _GxEScanR_makegxr(SEXP destSEXP, SEXP srcSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type dest(destSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type src(srcSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idx(idxSEXP);
    makegxr(dest, src, idx);
    return R_NilValue;
END_RCPP
}
// xrgwis2
void xrgwis2(arma::mat& xr1, arma::mat& xr2, arma::mat& xr3, arma::mat& xr4, arma::mat& xr5, const arma::uvec& idx, const arma::mat& src1, const arma::vec& src2);
RcppExport SEXP _GxEScanR_xrgwis2(SEXP xr1SEXP, SEXP xr2SEXP, SEXP xr3SEXP, SEXP xr4SEXP, SEXP xr5SEXP, SEXP idxSEXP, SEXP src1SEXP, SEXP src2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type xr1(xr1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr2(xr2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr3(xr3SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr4(xr4SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr5(xr5SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type src1(src1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type src2(src2SEXP);
    xrgwis2(xr1, xr2, xr3, xr4, xr5, idx, src1, src2);
    return R_NilValue;
END_RCPP
}
// lrtgwis2
void lrtgwis2(arma::vec& lrtg, arma::vec& lrtgxe, arma::vec& lrt2df, arma::vec& lrteg, arma::vec& lrt3df, arma::vec& lrtcase, arma::vec& lrtctrl, const arma::vec& loglike0, const arma::vec& loglhg, const arma::vec& loglhgxe, const arma::vec& loglheg, const arma::vec& loglhcase, const arma::vec& loglhctrl, const int offset);
RcppExport SEXP _GxEScanR_lrtgwis2(SEXP lrtgSEXP, SEXP lrtgxeSEXP, SEXP lrt2dfSEXP, SEXP lrtegSEXP, SEXP lrt3dfSEXP, SEXP lrtcaseSEXP, SEXP lrtctrlSEXP, SEXP loglike0SEXP, SEXP loglhgSEXP, SEXP loglhgxeSEXP, SEXP loglhegSEXP, SEXP loglhcaseSEXP, SEXP loglhctrlSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type lrtg(lrtgSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lrtgxe(lrtgxeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lrt2df(lrt2dfSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lrteg(lrtegSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lrt3df(lrt3dfSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lrtcase(lrtcaseSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lrtctrl(lrtctrlSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglike0(loglike0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglhg(loglhgSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglhgxe(loglhgxeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglheg(loglhegSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglhcase(loglhcaseSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglhctrl(loglhctrlSEXP);
    Rcpp::traits::input_parameter< const int >::type offset(offsetSEXP);
    lrtgwis2(lrtg, lrtgxe, lrt2df, lrteg, lrt3df, lrtcase, lrtctrl, loglike0, loglhg, loglhgxe, loglheg, loglhcase, loglhctrl, offset);
    return R_NilValue;
END_RCPP
}
// betagwis2
void betagwis2(arma::vec& betag, arma::vec& betagxe, arma::vec& betaeg, arma::vec& betacase, arma::vec& betactrl, const arma::mat& tmpbetag, const arma::mat& tmpbetagxe, const arma::mat& tmpbetaeg, const arma::mat& tmpbetacase, const arma::mat& tmpbetactrl, const double estddev, const int offset);
RcppExport SEXP _GxEScanR_betagwis2(SEXP betagSEXP, SEXP betagxeSEXP, SEXP betaegSEXP, SEXP betacaseSEXP, SEXP betactrlSEXP, SEXP tmpbetagSEXP, SEXP tmpbetagxeSEXP, SEXP tmpbetaegSEXP, SEXP tmpbetacaseSEXP, SEXP tmpbetactrlSEXP, SEXP estddevSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type betag(betagSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type betagxe(betagxeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type betaeg(betaegSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type betacase(betacaseSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type betactrl(betactrlSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tmpbetag(tmpbetagSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tmpbetagxe(tmpbetagxeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tmpbetaeg(tmpbetaegSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tmpbetacase(tmpbetacaseSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tmpbetactrl(tmpbetactrlSEXP);
    Rcpp::traits::input_parameter< const double >::type estddev(estddevSEXP);
    Rcpp::traits::input_parameter< const int >::type offset(offsetSEXP);
    betagwis2(betag, betagxe, betaeg, betacase, betactrl, tmpbetag, tmpbetagxe, tmpbetaeg, tmpbetacase, tmpbetactrl, estddev, offset);
    return R_NilValue;
END_RCPP
}
// makegxexr
void makegxexr(arma::mat& dest, const arma::mat& src1, const arma::mat& src2);
RcppExport SEXP _GxEScanR_makegxexr(SEXP destSEXP, SEXP src1SEXP, SEXP src2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type dest(destSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type src1(src1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type src2(src2SEXP);
    makegxexr(dest, src1, src2);
    return R_NilValue;
END_RCPP
}
// copybeta
void copybeta(arma::vec& dest, const arma::mat& src, int colnum, int startrow, int numrows);
RcppExport SEXP _GxEScanR_copybeta(SEXP destSEXP, SEXP srcSEXP, SEXP colnumSEXP, SEXP startrowSEXP, SEXP numrowsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type dest(destSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type src(srcSEXP);
    Rcpp::traits::input_parameter< int >::type colnum(colnumSEXP);
    Rcpp::traits::input_parameter< int >::type startrow(startrowSEXP);
    Rcpp::traits::input_parameter< int >::type numrows(numrowsSEXP);
    copybeta(dest, src, colnum, startrow, numrows);
    return R_NilValue;
END_RCPP
}
// calculatelrt
void calculatelrt(arma::vec& lrt, int idx1, int idx2, const arma::vec& loglike, double loglike0);
RcppExport SEXP _GxEScanR_calculatelrt(SEXP lrtSEXP, SEXP idx1SEXP, SEXP idx2SEXP, SEXP loglikeSEXP, SEXP loglike0SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type lrt(lrtSEXP);
    Rcpp::traits::input_parameter< int >::type idx1(idx1SEXP);
    Rcpp::traits::input_parameter< int >::type idx2(idx2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglike(loglikeSEXP);
    Rcpp::traits::input_parameter< double >::type loglike0(loglike0SEXP);
    calculatelrt(lrt, idx1, idx2, loglike, loglike0);
    return R_NilValue;
END_RCPP
}
// calculatelrtgxe
void calculatelrtgxe(arma::vec& lrtg, arma::vec& lrtgxe, arma::vec& lrt2df, int idx1, int idx2, const arma::vec& loglikeg, const arma::vec& loglikegxe, double loglike0);
RcppExport SEXP _GxEScanR_calculatelrtgxe(SEXP lrtgSEXP, SEXP lrtgxeSEXP, SEXP lrt2dfSEXP, SEXP idx1SEXP, SEXP idx2SEXP, SEXP loglikegSEXP, SEXP loglikegxeSEXP, SEXP loglike0SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type lrtg(lrtgSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lrtgxe(lrtgxeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lrt2df(lrt2dfSEXP);
    Rcpp::traits::input_parameter< int >::type idx1(idx1SEXP);
    Rcpp::traits::input_parameter< int >::type idx2(idx2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglikeg(loglikegSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type loglikegxe(loglikegxeSEXP);
    Rcpp::traits::input_parameter< double >::type loglike0(loglike0SEXP);
    calculatelrtgxe(lrtg, lrtgxe, lrt2df, idx1, idx2, loglikeg, loglikegxe, loglike0);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GxEScanR_initlslinreg", (DL_FUNC) &_GxEScanR_initlslinreg, 2},
    {"_GxEScanR_lslinreg", (DL_FUNC) &_GxEScanR_lslinreg, 14},
    {"_GxEScanR_initlslogreg", (DL_FUNC) &_GxEScanR_initlslogreg, 3},
    {"_GxEScanR_lslogreg", (DL_FUNC) &_GxEScanR_lslogreg, 18},
    {"_GxEScanR_readblock", (DL_FUNC) &_GxEScanR_readblock, 4},
    {"_GxEScanR_getdosages", (DL_FUNC) &_GxEScanR_getdosages, 7},
    {"_GxEScanR_stdmat", (DL_FUNC) &_GxEScanR_stdmat, 4},
    {"_GxEScanR_makegxr", (DL_FUNC) &_GxEScanR_makegxr, 3},
    {"_GxEScanR_xrgwis2", (DL_FUNC) &_GxEScanR_xrgwis2, 8},
    {"_GxEScanR_lrtgwis2", (DL_FUNC) &_GxEScanR_lrtgwis2, 14},
    {"_GxEScanR_betagwis2", (DL_FUNC) &_GxEScanR_betagwis2, 12},
    {"_GxEScanR_makegxexr", (DL_FUNC) &_GxEScanR_makegxexr, 3},
    {"_GxEScanR_copybeta", (DL_FUNC) &_GxEScanR_copybeta, 5},
    {"_GxEScanR_calculatelrt", (DL_FUNC) &_GxEScanR_calculatelrt, 5},
    {"_GxEScanR_calculatelrtgxe", (DL_FUNC) &_GxEScanR_calculatelrtgxe, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_GxEScanR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
