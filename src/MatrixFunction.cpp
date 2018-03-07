#include <cstdlib>
#include <cstring>
#include "MatrixFunctions.h"

// An important thing to remember here is that C and C++ are row dominant
//    n represents the number of rows
//    m represents the number of columns

void Transpose(const double *x, double *y, const unsigned int n, const unsigned int m) {
  unsigned int ui, uj;
  const double *x1, *x2;
  double *y1;
  
  x1 = x;
  y1 = y;
  for (ui = 0; ui < m; ++ui, ++x1) {
    x2 = x1;
    for (uj = 0; uj < n; ++uj, x2 += m, ++y1)
      *y1 = *x2;
  }
}

void XTransposeX(const double *x, double *y, const unsigned int n, const unsigned int m) {
  unsigned int ui, uj, uk;
  const double *xp1, *xp2, *xp3, *xp4;
  double *yp1, *yp2;
  
  xp1 = x;
  memset(y, 0, n*n*sizeof(double));
  yp1 = y;
  for (ui = 0; ui < n; ++ui, ++xp1) {
    yp1 += ui;
    yp2 = yp1;
    xp2 = xp1;
    for (uj = ui; uj < n; ++uj, ++xp2, ++yp1, yp2 += n) {
      xp3 = xp1;
      xp4 = xp2;
      for (uk = 0; uk < m; ++uk, xp3 += n, xp4 += n)
        *yp1 += *xp3 * * xp4;
      *yp2 = *yp1;
    }
  }
}

void XTransposeX(const double *x, double *y, const bool *_notMissing, const unsigned int n, const unsigned int m) {
  unsigned int ui, uj, uk;
  const double *xp1, *xp2, *xp3, *xp4;
  double *yp1, *yp2;
  const bool *notMissing;
  
  xp1 = x;
  memset(y, 0, n*n*sizeof(double));
  yp1 = y;
  for (ui = 0; ui < n; ++ui, ++xp1) {
    yp1 += ui;
    yp2 = yp1;
    xp2 = xp1;
    for (uj = ui; uj < n; ++uj, ++xp2, ++yp1, yp2 += n) {
      xp3 = xp1;
      xp4 = xp2;
      notMissing = _notMissing;
      for (uk = 0; uk < m; ++uk, xp3 += n, xp4 += n, ++notMissing) {
        if (*notMissing == true)
          *yp1 += *xp3 * * xp4;
      }
      *yp2 = *yp1;
    }
  }
}

void XTransposeX(const double *x, double *y, const bool *_notMissing, const unsigned int n, const unsigned int m, const unsigned int p) {
  unsigned int ui, uj, uk;
  const double *xp1, *xp2, *xp3, *xp4;
  double *yp1, *yp2;
  const bool *notMissing;
  
  xp1 = x;
  memset(y, 0, n*n*sizeof(double));
  yp1 = y;
  for (ui = 0; ui < n; ++ui, ++xp1) {
    yp1 += ui;
    yp2 = yp1;
    xp2 = xp1;
    for (uj = ui; uj < n; ++uj, ++xp2, ++yp1, yp2 += n) {
      xp3 = xp1;
      xp4 = xp2;
      notMissing = _notMissing;
      for (uk = 0; uk < m; ++uk, xp3 += p, xp4 += p, ++notMissing) {
        if (*notMissing == true)
          *yp1 += *xp3 * * xp4;
      }
      *yp2 = *yp1;
    }
  }
}

void XTransposeY(const double *x, const double *y, double *xty, const bool *_notMissing, const unsigned int n, const unsigned int m) {
  const double *xp1, *xp2, *yp;
  const bool *mp;
  double *xtyp;
  unsigned int ui, uj;
  
  xp1 = x;
  xtyp = xty;
  for (ui = 0; ui < n; ++ui, ++xp1, ++xtyp) {
    *xtyp = 0;
    xp2 = xp1;
    yp = y;
    mp = _notMissing;
    for (uj = 0; uj < m; ++uj, xp2 += n, ++yp, ++mp) {
      if (*mp == true)
        *xtyp += *xp2 * *yp;
    }
  }
}
// Cholesky upper-lower decomposition
// This only works for positive-definite matrices
// x is a nxn matrix, y will contain the upper lower decomposition
// The diagonal of y is for the upper triangular matrix,
// The lower triangular matrix has 1's along the diagonal.
int Cholesky(const double *x, double *y, int n, double evalue) {
  int i, j, k;
  double *cpi, *cpj;
  double *Lik, *Ljk;
  double *dpi, *dpk;
  double *yp;
  
  memmove((char *)y, (char *)x, n * n * sizeof(double));
  cpi = y;
  for (i = 0; i < n; ++i) {
    yp = cpi + i * n;
    dpk = y;
    Lik = cpi;
    for (k = 0; k < i; ++k) {
      *yp -= *Lik * *Lik * *dpk;
      Lik += n;
      dpk += n + 1;
    }
    // ??? Why is this check here?
    if (*dpk < evalue && *dpk != *dpk)
      return 1;
    //		if (*dpk < 1.0e-7)
    //			return 2;
    
    dpi = dpk;
    cpj = cpi;
    for (j = i + 1; j < n; ++j) {
      ++yp;
      ++cpj;
      Lik = cpi;
      Ljk = cpj;
      dpk = y;
      for (k = 0; k < i; ++k) {
        *yp -= *Lik * *Ljk * *dpk;
        Lik += n;
        Ljk += n;
        dpk += n + 1;
      }
      *yp /= *dpi;
    }
    ++cpi;
  }
  dpk = y;
  // Make sure diagonal values are valid
  for (i = 0; i < n; ++i) {
    if (*dpk < evalue || *dpk != *dpk)
      return 1;
    dpk += n + 1;
  }
  return 0;
}

void BackwardForward(double *x, double *y, int n) {
  int i, j, k;
  double *d1;
  double *xp1, *xp2;
  double *yp1, *yp2, *yp3;
  
  memset(y, 0, n * n * sizeof(double));
  yp1 = y;
  for (i = n; i; --i, yp1 += n + 1)
    *yp1 = 1.;
  yp1 = y;
  d1 = x;
  for (i = 0; i < n; ++i) {
    xp1 = d1;
    yp2 = yp1;
    for (j = i; j < n; ++j) {
      xp2 = xp1;
      yp3 = yp2;
      for (k = j; k; --k) {
        --yp3;
        xp2 -= n;
        *yp2 -= *xp2 * *yp3;
      }
      xp1 += n + 1;
      ++yp2;
    }
    yp1 += n + 1;
    d1 += n + 1;
  }
  d1 = x;
  yp1 = y;
  for (i = 0; i < n; ++i) {
    yp2 = yp1 + i;
    xp1 = d1;
    for (j = i; j < n; ++j) {
      *yp2 /= *xp1;
      ++yp2;
      xp1 += n + 1;
    }
    yp1 += n;
    d1 += n + 1;
  }
  
  d1 = x + n * n - 1;
  yp1 = y + n - 1;
  for (i = n; i; --i) {
    xp1 = d1;
    yp2 = yp1;
    for (j = 0; j < i; ++j) {
      xp2 = xp1;
      yp3 = yp1;
      for (k = j; k; --k) {
        *yp2 -= *xp2 * *yp3;
        --xp2;
        --yp3;
      }
      --yp2;
      xp1 -= n;
    }
    yp1 += n;
    //		d1 -= n + 1;
  }
  
}

int Invert(const double *x, double *ULDecomp, double *xinv, const unsigned int n, double evalue) {
  double *y1, *y2, *y3, *y4;
  unsigned int ui, uj;
  
  // Cholesky upper lower decomposition
  if (Cholesky(x, ULDecomp, n, evalue)) {
    return 1;
  }
  // Backward-forward substituion only fills in half the matrix
  BackwardForward(ULDecomp, xinv, n);
  // Fill in the second half of the matrix
  y1 = xinv;
  y2 = xinv;
  for (ui = 0; ui < n; ++ui, y1 += n, ++y2) {
    y3 = y1;
    y4 = y2;
    for (uj = 0; uj < ui; ++uj, ++y3, y4 += n)
      *y3 = *y4;
  }
  return 0;
}
