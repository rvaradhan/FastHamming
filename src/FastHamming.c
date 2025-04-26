#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>  // Required for OpenMP functions
#endif

// Fast popcount for 64-bit values
static inline int popcount64(uint64_t x) {
  return __builtin_popcountll(x);
}

// Compute pairwise Hamming distances, skipping diagonal since it's always zero
void compute_hamming_distances(const uint64_t* data, int n_rows, int n_words, int* result) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < n_rows; ++i) {
    // Set diagonal element to zero
    result[i + i * n_rows] = 0;
    for (int j = i + 1; j < n_rows; ++j) {
      int dist = 0;
      // XOR-packed words and accumulate popcounts
      const uint64_t* row_i = data + i * n_words;
      const uint64_t* row_j = data + j * n_words;
      for (int k = 0; k < n_words; ++k) {
        uint64_t xor_val = row_i[k] ^ row_j[k];
        dist += popcount64(xor_val);
      }
      // Fill symmetric entries
      result[i + j * n_rows] = dist;
      result[j + i * n_rows] = dist;
    }
  }
}

// R wrapper for .Call interface
SEXP c_hamming_distance(SEXP r_data,
                        SEXP r_n_rows,
                        SEXP r_n_cols,
                        SEXP r_nthreads) {
  // Retrieve dimensions
  int n_rows = INTEGER(r_n_rows)[0];
  int n_cols = INTEGER(r_n_cols)[0];
  int n_words = (n_cols + 63) / 64;

#ifdef _OPENMP
  // Thread setup: use user-specified or default to all available
  if (!isNull(r_nthreads)) {
    int user_threads = INTEGER(r_nthreads)[0];
    omp_set_dynamic(0);
    omp_set_num_threads(user_threads);
  } else {
    omp_set_num_threads(omp_get_max_threads());
  }
  // Report what OpenMP is actually doing
#pragma omp parallel
{
#pragma omp master
{
  int max_threads  = omp_get_max_threads();
  int num_procs    = omp_get_num_procs();
  int actual       = omp_get_num_threads();
  Rprintf("OpenMP reports: max_threads = %d, num_procs = %d, actual used = %d\n",
          max_threads, num_procs, actual);
}
}
#endif

// Allocate packed data buffer
uint64_t* packed_data = (uint64_t*)calloc((size_t)n_rows * n_words, sizeof(uint64_t));
if (!packed_data) {
  error("Memory allocation failed for packed data.");
}

// Pack binary matrix from INTEGER(r_data)
int* raw = INTEGER(r_data);
for (int i = 0; i < n_rows; ++i) {
  for (int j = 0; j < n_cols; ++j) {
    if (raw[i + j * n_rows]) {
      packed_data[i * n_words + (j >> 6)] |= ((uint64_t)1 << (j & 63));
    }
  }
}

// Allocate result matrix (INTSXP)
SEXP r_result = PROTECT(allocMatrix(INTSXP, n_rows, n_rows));
int* result = INTEGER(r_result);

// Compute distances
compute_hamming_distances(packed_data, n_rows, n_words, result);

// Cleanup and return
free(packed_data);
UNPROTECT(1);
return r_result;
}

