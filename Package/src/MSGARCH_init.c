// RegisteringDynamic Symbols
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls
*/

extern SEXP MSGARCH_getDelta(SEXP, SEXP);
extern SEXP MSGARCH_Viterbi(SEXP, SEXP, SEXP);
extern SEXP MSGARCH_EM_HMM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MSGARCH_EM_HMM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MSGARCH_EM_MM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MSGARCH_MapParameters_univ(SEXP, SEXP, SEXP);
extern SEXP MSGARCH_UnmapParameters_univ(SEXP, SEXP, SEXP);
extern SEXP MSGARCH_dUnivLike(SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP _rcpp_module_boot_MSgarch();
extern SEXP _rcpp_module_boot_tGARCH();
extern SEXP _rcpp_module_boot_sGARCH();
extern SEXP _rcpp_module_boot_dist();
extern SEXP _rcpp_module_boot_gjrGARCH();
extern SEXP _rcpp_module_boot_GAS();
extern SEXP _rcpp_module_boot_eGARCH();

static const R_CallMethodDef CallEntries[] = {
    {"MSGARCH_getDelta", (DL_FUNC) &MSGARCH_getDelta, 2},
    {"MSGARCH_Viterbi", (DL_FUNC) &MSGARCH_Viterbi, 3},
    {"MSGARCH_EM_HMM", (DL_FUNC) &MSGARCH_EM_HMM, 5},
    {"MSGARCH_EM_MM", (DL_FUNC) &MSGARCH_EM_MM, 5},
    {"MSGARCH_MapParameters_univ", (DL_FUNC) &MSGARCH_MapParameters_univ, 3},
    {"MSGARCH_UnmapParameters_univ", (DL_FUNC) &MSGARCH_UnmapParameters_univ, 3},
    {"MSGARCH_dUnivLike", (DL_FUNC) &MSGARCH_dUnivLike, 5},
	{"_rcpp_module_boot_MSgarch", (DL_FUNC) &_rcpp_module_boot_MSgarch,0},
	{"_rcpp_module_boot_tGARCH", (DL_FUNC) &_rcpp_module_boot_tGARCH,0},
	{"_rcpp_module_boot_sGARCH", (DL_FUNC) &_rcpp_module_boot_sGARCH,0},
	{"_rcpp_module_boot_dist", (DL_FUNC) &_rcpp_module_boot_dist,0},
	{"_rcpp_module_boot_gjrGARCH", (DL_FUNC) &_rcpp_module_boot_gjrGARCH,0},
	{"_rcpp_module_boot_GAS", (DL_FUNC) &_rcpp_module_boot_GAS,0},
	{"_rcpp_module_boot_eGARCH", (DL_FUNC) &_rcpp_module_boot_eGARCH,0},
    {NULL, NULL, 0}
};

void R_init_MSGARCH(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
