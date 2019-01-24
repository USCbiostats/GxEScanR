// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GetBinaryDosageInformation
Rcpp::List GetBinaryDosageInformation(const std::string& binaryDosageFilename, const unsigned int snSub, const unsigned int snSNPs);
RcppExport SEXP _GxEScanR_GetBinaryDosageInformation(SEXP binaryDosageFilenameSEXP, SEXP snSubSEXP, SEXP snSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type binaryDosageFilename(binaryDosageFilenameSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type snSub(snSubSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type snSNPs(snSNPsSEXP);
    rcpp_result_gen = Rcpp::wrap(GetBinaryDosageInformation(binaryDosageFilename, snSub, snSNPs));
    return rcpp_result_gen;
END_RCPP
}
// GxEScanC
int GxEScanC(Rcpp::List subjectData, Rcpp::List geneticInfo, std::string outputFilename, std::string skippedFilename, double minMaf, double geCutoff);
RcppExport SEXP _GxEScanR_GxEScanC(SEXP subjectDataSEXP, SEXP geneticInfoSEXP, SEXP outputFilenameSEXP, SEXP skippedFilenameSEXP, SEXP minMafSEXP, SEXP geCutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type subjectData(subjectDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type geneticInfo(geneticInfoSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilename(outputFilenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type skippedFilename(skippedFilenameSEXP);
    Rcpp::traits::input_parameter< double >::type minMaf(minMafSEXP);
    Rcpp::traits::input_parameter< double >::type geCutoff(geCutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(GxEScanC(subjectData, geneticInfo, outputFilename, skippedFilename, minMaf, geCutoff));
    return rcpp_result_gen;
END_RCPP
}
// GxEScanCSubset
int GxEScanCSubset(Rcpp::List subjectData, Rcpp::List geneticInfo, std::string outputFilename, std::string skippedFilename, double minMaf, double geCutoff, std::vector<int>& snpIndices);
RcppExport SEXP _GxEScanR_GxEScanCSubset(SEXP subjectDataSEXP, SEXP geneticInfoSEXP, SEXP outputFilenameSEXP, SEXP skippedFilenameSEXP, SEXP minMafSEXP, SEXP geCutoffSEXP, SEXP snpIndicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type subjectData(subjectDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type geneticInfo(geneticInfoSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilename(outputFilenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type skippedFilename(skippedFilenameSEXP);
    Rcpp::traits::input_parameter< double >::type minMaf(minMafSEXP);
    Rcpp::traits::input_parameter< double >::type geCutoff(geCutoffSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type snpIndices(snpIndicesSEXP);
    rcpp_result_gen = Rcpp::wrap(GxEScanCSubset(subjectData, geneticInfo, outputFilename, skippedFilename, minMaf, geCutoff, snpIndices));
    return rcpp_result_gen;
END_RCPP
}
// GxETest
int GxETest();
RcppExport SEXP _GxEScanR_GxETest() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(GxETest());
    return rcpp_result_gen;
END_RCPP
}
// GxEScanFreqC
int GxEScanFreqC(Rcpp::List subjectData, Rcpp::List geneticInfo, std::string outputFilename);
RcppExport SEXP _GxEScanR_GxEScanFreqC(SEXP subjectDataSEXP, SEXP geneticInfoSEXP, SEXP outputFilenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type subjectData(subjectDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type geneticInfo(geneticInfoSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilename(outputFilenameSEXP);
    rcpp_result_gen = Rcpp::wrap(GxEScanFreqC(subjectData, geneticInfo, outputFilename));
    return rcpp_result_gen;
END_RCPP
}
// GxEScanFreqC2
int GxEScanFreqC2(Rcpp::List subjectData, Rcpp::List geneticInfo, std::string outputFilename);
RcppExport SEXP _GxEScanR_GxEScanFreqC2(SEXP subjectDataSEXP, SEXP geneticInfoSEXP, SEXP outputFilenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type subjectData(subjectDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type geneticInfo(geneticInfoSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputFilename(outputFilenameSEXP);
    rcpp_result_gen = Rcpp::wrap(GxEScanFreqC2(subjectData, geneticInfo, outputFilename));
    return rcpp_result_gen;
END_RCPP
}
// Imp2toBDC
int Imp2toBDC(const Rcpp::List& imp2Info, const std::string& filename, int format, int subformat);
RcppExport SEXP _GxEScanR_Imp2toBDC(SEXP imp2InfoSEXP, SEXP filenameSEXP, SEXP formatSEXP, SEXP subformatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type imp2Info(imp2InfoSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type format(formatSEXP);
    Rcpp::traits::input_parameter< int >::type subformat(subformatSEXP);
    rcpp_result_gen = Rcpp::wrap(Imp2toBDC(imp2Info, filename, format, subformat));
    return rcpp_result_gen;
END_RCPP
}
// PlinkBinaryInfo
Rcpp::List PlinkBinaryInfo(std::string& geneticFile, std::string& mapFile, std::string& familyFile);
RcppExport SEXP _GxEScanR_PlinkBinaryInfo(SEXP geneticFileSEXP, SEXP mapFileSEXP, SEXP familyFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type geneticFile(geneticFileSEXP);
    Rcpp::traits::input_parameter< std::string& >::type mapFile(mapFileSEXP);
    Rcpp::traits::input_parameter< std::string& >::type familyFile(familyFileSEXP);
    rcpp_result_gen = Rcpp::wrap(PlinkBinaryInfo(geneticFile, mapFile, familyFile));
    return rcpp_result_gen;
END_RCPP
}
// OpenGxEOutFile
int OpenGxEOutFile(std::string& filename);
RcppExport SEXP _GxEScanR_OpenGxEOutFile(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(OpenGxEOutFile(filename));
    return rcpp_result_gen;
END_RCPP
}
// AppendGxEResults
int AppendGxEResults(std::string& filename, Rcpp::StringVector& snpID, Rcpp::StringVector& chromosome, Rcpp::IntegerVector& location, Rcpp::StringVector& refAllele, Rcpp::StringVector& altAllele, int numSub, int numCases, arma::mat& logLike, arma::mat& estimates, int length, double sigmaE);
RcppExport SEXP _GxEScanR_AppendGxEResults(SEXP filenameSEXP, SEXP snpIDSEXP, SEXP chromosomeSEXP, SEXP locationSEXP, SEXP refAlleleSEXP, SEXP altAlleleSEXP, SEXP numSubSEXP, SEXP numCasesSEXP, SEXP logLikeSEXP, SEXP estimatesSEXP, SEXP lengthSEXP, SEXP sigmaESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type snpID(snpIDSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type chromosome(chromosomeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type location(locationSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type refAllele(refAlleleSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type altAllele(altAlleleSEXP);
    Rcpp::traits::input_parameter< int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< int >::type numCases(numCasesSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type logLike(logLikeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type estimates(estimatesSEXP);
    Rcpp::traits::input_parameter< int >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaE(sigmaESEXP);
    rcpp_result_gen = Rcpp::wrap(AppendGxEResults(filename, snpID, chromosome, location, refAllele, altAllele, numSub, numCases, logLike, estimates, length, sigmaE));
    return rcpp_result_gen;
END_RCPP
}
// InitializeLRMod
int InitializeLRMod(int numRow, int numCol, arma::vec& y, arma::mat& xl, arma::vec& beta, arma::vec& score, arma::vec& w, arma::vec& wInv, arma::vec& yp, arma::vec& zt, arma::vec& k, arma::mat& ql, arma::mat& rtl, arma::vec& abx, arma::vec& expabx, arma::vec& expabxp1, arma::vec& expitabx, arma::vec& logLikelihood);
RcppExport SEXP _GxEScanR_InitializeLRMod(SEXP numRowSEXP, SEXP numColSEXP, SEXP ySEXP, SEXP xlSEXP, SEXP betaSEXP, SEXP scoreSEXP, SEXP wSEXP, SEXP wInvSEXP, SEXP ypSEXP, SEXP ztSEXP, SEXP kSEXP, SEXP qlSEXP, SEXP rtlSEXP, SEXP abxSEXP, SEXP expabxSEXP, SEXP expabxp1SEXP, SEXP expitabxSEXP, SEXP logLikelihoodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type numRow(numRowSEXP);
    Rcpp::traits::input_parameter< int >::type numCol(numColSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wInv(wInvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yp(ypSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zt(ztSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtl(rtlSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type abx(abxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expabx(expabxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expabxp1(expabxp1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expitabx(expitabxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type logLikelihood(logLikelihoodSEXP);
    rcpp_result_gen = Rcpp::wrap(InitializeLRMod(numRow, numCol, y, xl, beta, score, w, wInv, yp, zt, k, ql, rtl, abx, expabx, expabxp1, expitabx, logLikelihood));
    return rcpp_result_gen;
END_RCPP
}
// FitLRMod
int FitLRMod(int n, int p, arma::vec& y, arma::mat& xl, arma::mat& xr, arma::vec& beta0, arma::vec& score0, arma::vec& w, arma::vec& wInv, arma::vec& yp0, arma::vec& zt0, arma::vec& k0, arma::mat& ql, arma::mat& rtl, arma::vec& abx, arma::vec& expabx, arma::vec& expabxp1, arma::vec& expitabx, arma::vec& yp, arma::vec& zt, arma::vec& k, arma::vec& bt, arma::mat& xrw, arma::vec& beta, arma::vec& score, arma::vec& zb, arma::vec& bb, arma::mat& h, arma::mat& rtr, arma::mat& t, arma::mat& qr, arma::mat& rbr, arma::vec& logLikelihood);
RcppExport SEXP _GxEScanR_FitLRMod(SEXP nSEXP, SEXP pSEXP, SEXP ySEXP, SEXP xlSEXP, SEXP xrSEXP, SEXP beta0SEXP, SEXP score0SEXP, SEXP wSEXP, SEXP wInvSEXP, SEXP yp0SEXP, SEXP zt0SEXP, SEXP k0SEXP, SEXP qlSEXP, SEXP rtlSEXP, SEXP abxSEXP, SEXP expabxSEXP, SEXP expabxp1SEXP, SEXP expitabxSEXP, SEXP ypSEXP, SEXP ztSEXP, SEXP kSEXP, SEXP btSEXP, SEXP xrwSEXP, SEXP betaSEXP, SEXP scoreSEXP, SEXP zbSEXP, SEXP bbSEXP, SEXP hSEXP, SEXP rtrSEXP, SEXP tSEXP, SEXP qrSEXP, SEXP rbrSEXP, SEXP logLikelihoodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type score0(score0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wInv(wInvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yp0(yp0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zt0(zt0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtl(rtlSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type abx(abxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expabx(expabxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expabxp1(expabxp1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expitabx(expitabxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yp(ypSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zt(ztSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type bt(btSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xrw(xrwSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zb(zbSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtr(rtrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type qr(qrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rbr(rbrSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type logLikelihood(logLikelihoodSEXP);
    rcpp_result_gen = Rcpp::wrap(FitLRMod(n, p, y, xl, xr, beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl, abx, expabx, expabxp1, expitabx, yp, zt, k, bt, xrw, beta, score, zb, bb, h, rtr, t, qr, rbr, logLikelihood));
    return rcpp_result_gen;
END_RCPP
}
// ScanGenes
int ScanGenes(int n, int p, arma::vec& y, arma::mat& xl, arma::mat& xr, int numSNPs, double minMAF, arma::vec& beta0, arma::vec& score0, arma::vec& w, arma::vec& wInv, arma::vec& yp0, arma::vec& zt0, arma::vec& k0, arma::mat& ql, arma::mat& rtl, arma::vec& abx, arma::vec& expabx, arma::vec& expabxp1, arma::vec& expitabx, arma::vec& yp, arma::vec& zt, arma::vec& k, arma::vec& bt, arma::mat& xrw1, arma::vec& beta1, arma::vec& score1, arma::vec& zb1, arma::vec& bb1, arma::mat& h1, arma::mat& rtr1, arma::mat& t1, arma::mat& qr1, arma::mat& rbr1, arma::vec& logLikelihood1, arma::mat& xrw2, arma::vec& beta2, arma::vec& score2, arma::vec& zb2, arma::vec& bb2, arma::mat& h2, arma::mat& rtr2, arma::mat& t2, arma::mat& qr2, arma::mat& rbr2, arma::vec& logLikelihood2, arma::mat& xr1, arma::mat& xr2, arma::vec logLikelihood0, arma::mat& logLikelihoods, arma::mat& estimates);
RcppExport SEXP _GxEScanR_ScanGenes(SEXP nSEXP, SEXP pSEXP, SEXP ySEXP, SEXP xlSEXP, SEXP xrSEXP, SEXP numSNPsSEXP, SEXP minMAFSEXP, SEXP beta0SEXP, SEXP score0SEXP, SEXP wSEXP, SEXP wInvSEXP, SEXP yp0SEXP, SEXP zt0SEXP, SEXP k0SEXP, SEXP qlSEXP, SEXP rtlSEXP, SEXP abxSEXP, SEXP expabxSEXP, SEXP expabxp1SEXP, SEXP expitabxSEXP, SEXP ypSEXP, SEXP ztSEXP, SEXP kSEXP, SEXP btSEXP, SEXP xrw1SEXP, SEXP beta1SEXP, SEXP score1SEXP, SEXP zb1SEXP, SEXP bb1SEXP, SEXP h1SEXP, SEXP rtr1SEXP, SEXP t1SEXP, SEXP qr1SEXP, SEXP rbr1SEXP, SEXP logLikelihood1SEXP, SEXP xrw2SEXP, SEXP beta2SEXP, SEXP score2SEXP, SEXP zb2SEXP, SEXP bb2SEXP, SEXP h2SEXP, SEXP rtr2SEXP, SEXP t2SEXP, SEXP qr2SEXP, SEXP rbr2SEXP, SEXP logLikelihood2SEXP, SEXP xr1SEXP, SEXP xr2SEXP, SEXP logLikelihood0SEXP, SEXP logLikelihoodsSEXP, SEXP estimatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPs(numSNPsSEXP);
    Rcpp::traits::input_parameter< double >::type minMAF(minMAFSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type score0(score0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type wInv(wInvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yp0(yp0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zt0(zt0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtl(rtlSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type abx(abxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expabx(expabxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expabxp1(expabxp1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type expitabx(expitabxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yp(ypSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zt(ztSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type bt(btSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xrw1(xrw1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type score1(score1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zb1(zb1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type bb1(bb1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtr1(rtr1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type qr1(qr1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rbr1(rbr1SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type logLikelihood1(logLikelihood1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xrw2(xrw2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type score2(score2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zb2(zb2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type bb2(bb2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtr2(rtr2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type qr2(qr2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rbr2(rbr2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type logLikelihood2(logLikelihood2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr1(xr1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr2(xr2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type logLikelihood0(logLikelihood0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type logLikelihoods(logLikelihoodsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type estimates(estimatesSEXP);
    rcpp_result_gen = Rcpp::wrap(ScanGenes(n, p, y, xl, xr, numSNPs, minMAF, beta0, score0, w, wInv, yp0, zt0, k0, ql, rtl, abx, expabx, expabxp1, expitabx, yp, zt, k, bt, xrw1, beta1, score1, zb1, bb1, h1, rtr1, t1, qr1, rbr1, logLikelihood1, xrw2, beta2, score2, zb2, bb2, h2, rtr2, t2, qr2, rbr2, logLikelihood2, xr1, xr2, logLikelihood0, logLikelihoods, estimates));
    return rcpp_result_gen;
END_RCPP
}
// GetLocations
int GetLocations(Rcpp::IntegerVector& x, Rcpp::NumericVector& y, double fileSize, Rcpp::IntegerVector& bufferSize);
RcppExport SEXP _GxEScanR_GetLocations(SEXP xSEXP, SEXP ySEXP, SEXP fileSizeSEXP, SEXP bufferSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type fileSize(fileSizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type bufferSize(bufferSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLocations(x, y, fileSize, bufferSize));
    return rcpp_result_gen;
END_RCPP
}
// WriteLocations
int WriteLocations(Rcpp::NumericVector& y, Rcpp::IntegerVector& x);
RcppExport SEXP _GxEScanR_WriteLocations(SEXP ySEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteLocations(y, x));
    return rcpp_result_gen;
END_RCPP
}
// GetSections
int GetSections(Rcpp::NumericVector& y, Rcpp::IntegerVector& snpSection, Rcpp::NumericVector& fileLocation, Rcpp::IntegerVector& snpLocation, Rcpp::IntegerVector& bufferSize);
RcppExport SEXP _GxEScanR_GetSections(SEXP ySEXP, SEXP snpSectionSEXP, SEXP fileLocationSEXP, SEXP snpLocationSEXP, SEXP bufferSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type snpSection(snpSectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type fileLocation(fileLocationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type snpLocation(snpLocationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type bufferSize(bufferSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSections(y, snpSection, fileLocation, snpLocation, bufferSize));
    return rcpp_result_gen;
END_RCPP
}
// ReadSNP
int ReadSNP(Rcpp::IntegerVector& snpNumber, Rcpp::IntegerVector& subjectNumber, std::string& filename, Rcpp::IntegerVector& format, int numSub, int numSNPs, int bufferSize, Rcpp::IntegerVector& buffer, int sections, Rcpp::IntegerVector& snpSection, Rcpp::NumericVector& fileLocation, Rcpp::IntegerVector& snpLocation, Rcpp::IntegerVector& currentSection, Rcpp::NumericVector& dosage, Rcpp::NumericVector& p0, Rcpp::NumericVector& p1, Rcpp::NumericVector& p2, Rcpp::NumericMatrix& values);
RcppExport SEXP _GxEScanR_ReadSNP(SEXP snpNumberSEXP, SEXP subjectNumberSEXP, SEXP filenameSEXP, SEXP formatSEXP, SEXP numSubSEXP, SEXP numSNPsSEXP, SEXP bufferSizeSEXP, SEXP bufferSEXP, SEXP sectionsSEXP, SEXP snpSectionSEXP, SEXP fileLocationSEXP, SEXP snpLocationSEXP, SEXP currentSectionSEXP, SEXP dosageSEXP, SEXP p0SEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type snpNumber(snpNumberSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type subjectNumber(subjectNumberSEXP);
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type format(formatSEXP);
    Rcpp::traits::input_parameter< int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPs(numSNPsSEXP);
    Rcpp::traits::input_parameter< int >::type bufferSize(bufferSizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type buffer(bufferSEXP);
    Rcpp::traits::input_parameter< int >::type sections(sectionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type snpSection(snpSectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type fileLocation(fileLocationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type snpLocation(snpLocationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type currentSection(currentSectionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadSNP(snpNumber, subjectNumber, filename, format, numSub, numSNPs, bufferSize, buffer, sections, snpSection, fileLocation, snpLocation, currentSection, dosage, p0, p1, p2, values));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GxEScanR_GetBinaryDosageInformation", (DL_FUNC) &_GxEScanR_GetBinaryDosageInformation, 3},
    {"_GxEScanR_GxEScanC", (DL_FUNC) &_GxEScanR_GxEScanC, 6},
    {"_GxEScanR_GxEScanCSubset", (DL_FUNC) &_GxEScanR_GxEScanCSubset, 7},
    {"_GxEScanR_GxETest", (DL_FUNC) &_GxEScanR_GxETest, 0},
    {"_GxEScanR_GxEScanFreqC", (DL_FUNC) &_GxEScanR_GxEScanFreqC, 3},
    {"_GxEScanR_GxEScanFreqC2", (DL_FUNC) &_GxEScanR_GxEScanFreqC2, 3},
    {"_GxEScanR_Imp2toBDC", (DL_FUNC) &_GxEScanR_Imp2toBDC, 4},
    {"_GxEScanR_PlinkBinaryInfo", (DL_FUNC) &_GxEScanR_PlinkBinaryInfo, 3},
    {"_GxEScanR_OpenGxEOutFile", (DL_FUNC) &_GxEScanR_OpenGxEOutFile, 1},
    {"_GxEScanR_AppendGxEResults", (DL_FUNC) &_GxEScanR_AppendGxEResults, 12},
    {"_GxEScanR_InitializeLRMod", (DL_FUNC) &_GxEScanR_InitializeLRMod, 18},
    {"_GxEScanR_FitLRMod", (DL_FUNC) &_GxEScanR_FitLRMod, 33},
    {"_GxEScanR_ScanGenes", (DL_FUNC) &_GxEScanR_ScanGenes, 51},
    {"_GxEScanR_GetLocations", (DL_FUNC) &_GxEScanR_GetLocations, 4},
    {"_GxEScanR_WriteLocations", (DL_FUNC) &_GxEScanR_WriteLocations, 2},
    {"_GxEScanR_GetSections", (DL_FUNC) &_GxEScanR_GetSections, 5},
    {"_GxEScanR_ReadSNP", (DL_FUNC) &_GxEScanR_ReadSNP, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_GxEScanR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
