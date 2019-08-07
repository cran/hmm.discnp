#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(getgl)(double *, int *, int *,double *, double *, double *,int *, int *, int *, double *, int *, double *, double *, double *, double *, double *, double *);
extern void F77_NAME(gethgl)(double *, int *, int *, double *, double *, double *, double *, int *, int *, int *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
extern void F77_NAME(getl)(double *, double *, double *, int *, int *, double *, double *, double *);
extern void F77_NAME(recurse)(double *, double *, double *, int *, double *, int *, int *, int *, int *, double *, double *, int *, int *, double *, double *, double *, double *, double *);

static const R_FortranMethodDef FortranEntries[] = {
    {"getgl",   (DL_FUNC) &F77_NAME(getgl),   17},
    {"gethgl",  (DL_FUNC) &F77_NAME(gethgl),  25},
    {"getl",    (DL_FUNC) &F77_NAME(getl),     8},
    {"recurse", (DL_FUNC) &F77_NAME(recurse), 18},
    {NULL, NULL, 0}
};

void R_init_hmm_discnp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
