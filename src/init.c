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
extern SEXP _support_energycrit(SEXP, SEXP);
extern SEXP _support_gamma_eval(SEXP, SEXP);
extern SEXP _support_grad_qsp(SEXP, SEXP, SEXP);
extern SEXP _support_obj_qsp(SEXP, SEXP, SEXP);
extern SEXP _support_omega(SEXP, SEXP, SEXP);
extern SEXP _support_opt_hess(SEXP, SEXP);
extern SEXP _support_printBar(SEXP);
extern SEXP _support_psp_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _support_psp_mi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _support_psp_seq_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _support_sp_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _support_sp_seq_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _support_starL2cpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_support_csample_num",     (DL_FUNC) &_support_csample_num,      3},
    {"_support_energy_norm_cpp", (DL_FUNC) &_support_energy_norm_cpp,  3},
    {"_support_energycrit",      (DL_FUNC) &_support_energycrit,       2},
    {"_support_gamma_eval",      (DL_FUNC) &_support_gamma_eval,       2},
    {"_support_grad_qsp",        (DL_FUNC) &_support_grad_qsp,         3},
    {"_support_obj_qsp",         (DL_FUNC) &_support_obj_qsp,          3},
    {"_support_omega",           (DL_FUNC) &_support_omega,            3},
    {"_support_opt_hess",        (DL_FUNC) &_support_opt_hess,         2},
    {"_support_printBar",        (DL_FUNC) &_support_printBar,         1},
    {"_support_psp_cpp",         (DL_FUNC) &_support_psp_cpp,         16},
    {"_support_psp_mi",          (DL_FUNC) &_support_psp_mi,           6},
    {"_support_psp_seq_cpp",     (DL_FUNC) &_support_psp_seq_cpp,     18},
    {"_support_sp_cpp",          (DL_FUNC) &_support_sp_cpp,          15},
    {"_support_sp_seq_cpp",      (DL_FUNC) &_support_sp_seq_cpp,      13},
    {"_support_starL2cpp",       (DL_FUNC) &_support_starL2cpp,        3},
    {NULL, NULL, 0}
};

void R_init_support(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
