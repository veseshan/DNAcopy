#include <R.h>
void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { GetRNGstate(); }
double F77_SUB(dunif)(void) { return unif_rand(); }
