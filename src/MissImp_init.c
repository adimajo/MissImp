#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void classForest(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void regForest(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"classForest", (DL_FUNC) &classForest, 24},
  {"regForest",   (DL_FUNC) &regForest,   21},
  {NULL, NULL, 0}
};

void R_init_MissImp(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}