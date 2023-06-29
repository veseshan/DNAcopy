#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  NOTE: Fortran is not case sensitive; But in this file F77_NAME
        has to be in lower case to prevent "undefined symbol: smoothLR_" error
*/

/* .Fortran calls */
extern void F77_NAME(bsegci)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bsegp)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(esegp)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(fndcpt)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(getbdry)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(prune)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(smoothlr)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(wfindcpt)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"bsegci",   (DL_FUNC) &F77_NAME(bsegci),    9},
    {"bsegp",    (DL_FUNC) &F77_NAME(bsegp),     6},
    {"esegp",    (DL_FUNC) &F77_NAME(esegp),     7},
    {"fndcpt",   (DL_FUNC) &F77_NAME(fndcpt),   18},
    {"getbdry",  (DL_FUNC) &F77_NAME(getbdry),   7},
    {"prune",    (DL_FUNC) &F77_NAME(prune),    10},
    {"smoothLR", (DL_FUNC) &F77_NAME(smoothlr),  8},
    {"wfindcpt", (DL_FUNC) &F77_NAME(fndcpt),   21},
    {NULL, NULL, 0}
};

void R_init_clinfun(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
