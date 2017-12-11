#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _support_csample_num(SEXP, SEXP, SEXP);
extern SEXP _support_energy_norm_cpp(SEXP, SEXP, SEXP);
extern SEXP _support_energycrit(SEXP, SEXP, SEXP, SEXP);
extern SEXP _support_grad_qsp(SEXP, SEXP, SEXP);
extern SEXP _support_obj_qsp(SEXP, SEXP, SEXP);
extern SEXP _support_sp_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _support_starL2cpp(SEXP, SEXP, SEXP);
extern SEXP _support_thincpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_support_csample_num",     (DL_FUNC) &_support_csample_num,      3},
    {"_support_energy_norm_cpp", (DL_FUNC) &_support_energy_norm_cpp,  3},
    {"_support_energycrit",      (DL_FUNC) &_support_energycrit,       4},
    {"_support_grad_qsp",        (DL_FUNC) &_support_grad_qsp,         3},
    {"_support_obj_qsp",         (DL_FUNC) &_support_obj_qsp,          3},
    {"_support_sp_cpp",          (DL_FUNC) &_support_sp_cpp,          12},
    {"_support_starL2cpp",       (DL_FUNC) &_support_starL2cpp,        3},
    {"_support_thincpp",         (DL_FUNC) &_support_thincpp,         11},
    {NULL, NULL, 0}
};

void R_init_support(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
