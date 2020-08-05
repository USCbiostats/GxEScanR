// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

extern const int NUMBEROFBASES = 3;
// 0x7ffe is 32,767 or 2^16 - 1
// 0xfffe is 65,534 or 2^32 - 1
// 0x2710 is 10,000
extern const unsigned short USBASE[NUMBEROFBASES] = {
  0x7ffe, // Used for format 1.1
  0xfffe, // Used for format 1.2
  0x2710  // Used for all other formats (0x2710 = 10000)
};

// Values the short integers are multiplied by to get dosage and genetic
// probabilities
extern const double DBASE[NUMBEROFBASES] = {
  1. / USBASE[0],
  1. / USBASE[1],
  1. / USBASE[2]
};

// Read the block
// [[Rcpp::export]]
void readblock(Rcpp::StringVector &filename,
               const double blkloc,
               const double blkbytes,
               arma::ivec &blkbuffer) {
  std::ifstream bdfile;
  
  bdfile.open(filename[0], std::ios_base::in | std::ios_base::binary);
  if (!bdfile.good()) {
    Rcpp::Rcerr << "Failed to open file" << std::endl;
    return;
  }
  
  bdfile.seekg(blkloc);
  bdfile.read((char *)blkbuffer.memptr(), blkbytes);
  if (bdfile.fail())
    Rcpp::Rcerr << "Read error" << std::endl;
  
  bdfile.close();
}

// [[Rcpp::export]]
void getdosages(arma::mat &dosages,
                arma::ivec &blkbuffer,
                double fileloc,
                arma::vec &indices,
                int firstsnp,
                int lastsnp,
                int base) {
  double mult;
  unsigned short usmask;
  long long fileoffset;
  char *chptr;

  int currentsnp;
  double *snpoffset;
  double *d;
  unsigned short *usptr;
  unsigned int ui;
  
  mult = DBASE[base];
  usmask = (base == 0) ? 0xffff : 0x7fff;
  chptr = (char *)blkbuffer.memptr();
  fileoffset = (long long)fileloc;
  
  currentsnp = firstsnp - 1;
  snpoffset = indices.memptr() + currentsnp;
  d = dosages.memptr();
  for (; currentsnp < lastsnp; ++currentsnp, ++snpoffset) {
    usptr = (unsigned short *)(chptr + (long long)*snpoffset - fileoffset);
    for (ui = 0; ui < dosages.n_rows; ++ui, ++usptr, ++d)
      *d = (*usptr & usmask) * mult;
  }
}

// [[Rcpp::export]]
void stdmat(arma::mat &m1,
            arma::mat &m2,
            arma::rowvec &means,
            arma::rowvec &stddevs) {
  unsigned int ui;
  
  m2 = m1;
  means = mean(m2, 0);
  stddevs = stddev(m2, 0, 0);
  for (ui = 1; ui < m2.n_cols; ++ui) {
    m2.col(ui) -= means[ui];
    m2.col(ui) /= stddevs[ui];
  }
}

// [[Rcpp::export]]
void makegxr(arma::mat &dest,
             const arma::mat &src,
             const arma::uvec &idx) {
  int ncol;
  int nrowdest, nrowsrc;
  double *d1;
  const double *d2;
  const unsigned int *r;
  int ui, uj;

  ncol = dest.n_cols;
  nrowdest = dest.n_rows;
  nrowsrc = src.n_rows;
  d1 = dest.memptr();
  d2 = src.memptr();
  for (ui = 0; ui < ncol; ++ui, d2 += nrowsrc) {
    r = idx.memptr();
    for (uj = 0; uj < nrowdest; ++uj, ++r, ++d1)
      *d1 = *(d2 + *r);
  }
}

// [[Rcpp::export]]
void xrgwis2(arma::mat &xr1,
             arma::mat &xr2,
             arma::mat &xr3,
             arma::mat &xr4,
             arma::mat &xr5,
             const arma::uvec &idx,
             const arma::mat &src1,
             const arma::vec &src2) {
  double *d1, *d2a, *d2b, *d3, *d4, *d5;
  const double *cd1, *cd2;
  const unsigned int *idxp;
  double dose;
  unsigned int ncol;
  unsigned int nrowdest, nrowsrc;
  unsigned int ui, uj;
  
  ncol = src1.n_cols;
  nrowsrc = src1.n_rows;
  nrowdest = src2.n_elem;
  
  d1 = xr1.memptr();
  d2a = xr2.memptr();
  d2b = xr2.memptr() + nrowsrc;
  d3 = xr3.memptr();
  d4 = xr4.memptr();
  d5 = xr5.memptr();
  cd1 = src1.memptr();
  for (ui = 0 ; ui < ncol; ++ui, d2a += nrowdest, d2b += nrowdest, cd1 += nrowsrc) {
    cd2 = src2.memptr();
    idxp = idx.memptr();
    for (uj = 0; uj < nrowdest; ++uj, ++d1, ++d2a, ++d2b, ++d3, ++cd2, ++idxp) {
      dose = *(cd1 + *idxp);
//      Rcpp::Rcout << *(cd1 + *idxp) << '\t' << *idxp << '\t' << dose << std::endl;
      *d1 = dose;
      *d2a = dose;
      *d2b = dose * *cd2;
      *d3 = dose;
      if (uj < xr5.n_rows) {
        *d5 = dose;
        ++d5;
      } else {
        *d4 = dose;
        ++d4;
      }
    }
  }
//  Rcpp::Rcout << std::endl;
}

// [[Rcpp::export]]
void lrtgwis2(arma::vec &lrtg,
              arma::vec &lrtgxe,
              arma::vec &lrt2df,
              arma::vec &lrteg,
              arma::vec &lrt3df,
              arma::vec &lrtcase,
              arma::vec &lrtctrl,
              const arma::vec &loglike0,
              const arma::vec &loglhg,
              const arma::vec &loglhgxe,
              const arma::vec &loglheg,
              const arma::vec &loglhcase,
              const arma::vec &loglhctrl,
              const int offset) {
  double *d1, *d2, *d3, *d4, *d5, *d6, *d7;
  const double *cd1, *cd2, *cd3, *cd4, *cd5;
  unsigned int nelem;
  unsigned int ui;

  nelem = loglhg.n_elem;  
  d1 = lrtg.memptr();
  d2 = lrtgxe.memptr();
  d3 = lrt2df.memptr();
  d4 = lrteg.memptr();
  d5 = lrt3df.memptr();
  d6 = lrtcase.memptr();
  d7 = lrtctrl.memptr();
  cd1 = loglhg.memptr() + (offset - 1);
  cd2 = loglhgxe.memptr() + (offset - 1);
  cd3 = loglheg.memptr() + (offset - 1);
  cd4 = loglhcase.memptr() + (offset - 1);
  cd5 = loglhctrl.memptr() + (offset - 1);
  for (ui = 0; ui < nelem; ++ui, ++d1, ++d2, ++d3, ++d4, ++d5, ++d6, ++d7,
       ++cd1, ++cd2, ++cd3, ++cd4, ++cd5) {
    *d1 = 2 * (*cd1 - loglike0(0));
    *d3 = 2 * (*cd2 - loglike0(1));
    *d2 = *d3 - *d1;
    *d4 = 2 *(*cd3 - loglike0(2));
    *d5 = *d4 + *d3;
    *d6 = 2 * (*cd4 - loglike0(3));
    *d7 = 2 * (*cd5 - loglike0(4));
  }
}

// [[Rcpp::export]]
void betagwis2(arma::vec &betag,
               arma::vec &betagxe,
               arma::vec &betaeg,
               arma::vec &betacase,
               arma::vec &betactrl,
               const arma::mat &tmpbetag,
               const arma::mat &tmpbetagxe,
               const arma::mat &tmpbetaeg,
               const arma::mat &tmpbetacase,
               const arma::mat &tmpbetactrl,
               const double estddev,
               const int offset) {
  double *d1, *d2, *d3, *d4, *d5;
  const double *cd1, *cd2, *cd3, *cd4, *cd5;
  unsigned int nelem;
  unsigned int ui;
  
  nelem = tmpbetag.n_rows;
  d1 = betag.memptr();
  d2 = betagxe.memptr();
  d3 = betaeg.memptr();
  d4 = betacase.memptr();
  d5 = betactrl.memptr();
  cd1 = tmpbetag.memptr() + (offset - 1);
  cd2 = tmpbetagxe.memptr() + nelem + (offset - 1);
  cd3 = tmpbetaeg.memptr() + (offset - 1);
  cd4 = tmpbetacase.memptr() + (offset - 1);
  cd5 = tmpbetactrl.memptr() + (offset - 1);
  for (ui = 0; ui < nelem; ++ui, ++d1, ++d2, ++d3, ++d4, ++d5,
       ++cd1, ++cd2, ++cd3, ++cd4, ++cd5) {
    *d1 = *cd1;
    *d2 = *cd2 / estddev;
    *d3 = *cd3;
    *d4 = *cd4;
    *d5 = *cd5;
  }
}
// [[Rcpp::export]]
void makegxexr(arma::mat &dest,
               const arma::mat &src1,
               const arma::mat &src2) {
  int ncol;
  int nrows;
  double *d1, *d2;
  const double *d3a, *d3b, *d4;
  int ui, uj;
  
  ncol = src1.n_cols;
  nrows = src1.n_rows;
  d1 = dest.memptr();
  d2 = d1 + nrows;
  d3a = src2.memptr() + (src2.n_cols - 1) * src2.n_rows;
  d4 = src1.memptr();
  for (ui = 0; ui < ncol; ++ui, d1 += nrows, d2 += nrows) {
    d3b = d3a;
    for (uj = 0; uj < nrows; ++uj, ++d1, ++d2, ++d3b, ++d4) {
      *d1 = *d4;
      *d2 = *d3b * *d4;
    }
  }
}

// [[Rcpp::export]]
void copybeta(arma::vec &dest,
             const arma::mat &src,
             int colnum,
             int startrow,
             int numrows) {
  double *d;
  const double *s;
  
  d = dest.memptr() + startrow - 1;
  s = src.memptr() + (colnum - 1) * src.n_rows;
  memmove(d, s, numrows*sizeof(double));
}

// [[Rcpp::export]]
void copysubmat(arma::mat &dest,
                int drow1,
                int drow2,
                int dcol1,
                int dcol2,
                const arma::mat &src,
                int srow1,
                int srow2,
                int scol1,
                int scol2) {
  dest.submat(drow1 - 1, dcol1 - 1, drow2 - 1, dcol2 - 1) =
    src.submat(srow1 - 1, scol1 - 1, srow2 - 1, scol2 - 1);
}

// [[Rcpp::export]]
void calculatelrt(arma::vec &lrt,
                  int idx1,
                  int idx2,
                  const arma::vec &loglike,
                  double loglike0) {
  lrt.subvec(idx1-1, idx2-1) = 2 * (loglike.subvec(0, idx2-idx1) - loglike0);
}

// [[Rcpp::export]]
void calculatelrtgxe(arma::vec &lrtg,
                     arma::vec &lrtgxe,
                     arma::vec &lrt2df,
                    int idx1,
                    int idx2,
                    const arma::vec &loglikeg,
                    const arma::vec &loglikegxe,
                    double loglike0) {
  lrtg.subvec(idx1 - 1, idx2 - 1) =
    2 * (loglikeg.subvec(0, idx2 - idx1) - loglike0);
  lrt2df.subvec(idx1-1, idx2-1) =
    2 * (loglikegxe.subvec(0, idx2 - idx1) - loglike0);
  lrtgxe.subvec(idx1 - 1, idx2 - 1) =
    lrt2df.subvec(idx1 - 1, idx2 - 1) - lrtg.subvec(idx1 - 1, idx2 - 1);
}