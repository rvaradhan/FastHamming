#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare with 3 arguments
extern SEXP c_hamming_distance(SEXP r_data, SEXP r_n_rows, SEXP r_n_cols, SEXP r_nthreads);

// Register with correct argument count
static const R_CallMethodDef CallEntries[] = {
  {"c_hamming_distance", (DL_FUNC) &c_hamming_distance, 4},
  {NULL, NULL, 0}
};

void R_init_FastHamming(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
