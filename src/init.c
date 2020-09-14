#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _GxEScanR_betagwis2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_calculatelrt(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_calculatelrtgxe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_copybeta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_copysubmat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_getdosages(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_initlslinreg(SEXP, SEXP);
extern SEXP _GxEScanR_initlslogreg(SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_lrtgwis2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_lslinreg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_lslogreg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_makegxexr(SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_makegxr(SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_readblock(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_stdmat(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GxEScanR_xrgwis2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_GxEScanR_betagwis2",       (DL_FUNC) &_GxEScanR_betagwis2,       12},
  {"_GxEScanR_calculatelrt",    (DL_FUNC) &_GxEScanR_calculatelrt,     5},
  {"_GxEScanR_calculatelrtgxe", (DL_FUNC) &_GxEScanR_calculatelrtgxe,  8},
  {"_GxEScanR_copybeta",        (DL_FUNC) &_GxEScanR_copybeta,         5},
  {"_GxEScanR_copysubmat",      (DL_FUNC) &_GxEScanR_copysubmat,      10},
  {"_GxEScanR_getdosages",      (DL_FUNC) &_GxEScanR_getdosages,       7},
  {"_GxEScanR_initlslinreg",    (DL_FUNC) &_GxEScanR_initlslinreg,     2},
  {"_GxEScanR_initlslogreg",    (DL_FUNC) &_GxEScanR_initlslogreg,     3},
  {"_GxEScanR_lrtgwis2",        (DL_FUNC) &_GxEScanR_lrtgwis2,        14},
  {"_GxEScanR_lslinreg",        (DL_FUNC) &_GxEScanR_lslinreg,        14},
  {"_GxEScanR_lslogreg",        (DL_FUNC) &_GxEScanR_lslogreg,        18},
  {"_GxEScanR_makegxexr",       (DL_FUNC) &_GxEScanR_makegxexr,        3},
  {"_GxEScanR_makegxr",         (DL_FUNC) &_GxEScanR_makegxr,          3},
  {"_GxEScanR_readblock",       (DL_FUNC) &_GxEScanR_readblock,        4},
  {"_GxEScanR_stdmat",          (DL_FUNC) &_GxEScanR_stdmat,           4},
  {"_GxEScanR_xrgwis2",         (DL_FUNC) &_GxEScanR_xrgwis2,          8},
  {NULL, NULL, 0}
};

void R_init_GxEScanR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}