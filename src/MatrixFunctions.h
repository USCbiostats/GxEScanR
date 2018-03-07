#ifndef MATRIXFUNCTIONS_H
#define MATRIXFUNCTIONS_H 1

void Transpose(const double *x, double *y, const unsigned int n, const unsigned int m);
void XTransposeX(const double *x, double *y, const unsigned int n, const unsigned int m);
void XTransposeX(const double *x, double *y, const bool *_notMissing, const unsigned int n, const unsigned int m);
void XTransposeX(const double *x, double *y, const bool *_notMissing, const unsigned int n, const unsigned int m, const unsigned int p);
void XTransposeY(const double *x, const double *y, double *xty, const bool *_notMissing, const unsigned int n, const unsigned int m);
int Cholesky(const double *x, double *y, int n, double evalue = 1.0e-15);
void BackwardForward(double *x, double *y, int n);
int Invert(const double *x, double *ULDecomp, double *xinv, const unsigned int n, double evalue = 1.0e-15);

#endif
