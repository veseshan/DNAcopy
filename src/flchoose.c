#include <R.h>
#include <Rmath.h>
/* Fortran function for log-combinations */
double F77_SUB(flchoose)(double *n, double *k) { return lchoose(*n, *k); }
