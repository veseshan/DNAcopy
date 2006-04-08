#include <R.h>
#include <Rmath.h>
/* Fortran function for hypergeometric CDF */
double F77_SUB(fphypr)(double *i, double *m, double *n, double *k) { return phyper(*i, *m, *n, *k, 1, 0); }
