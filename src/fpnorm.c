#include <R.h>
#include <Rmath.h>
/* Fortran function for standard normal CDF */
double F77_SUB(fpnorm)(double *x) { return pnorm(*x, 0, 1, 1, 0); }
