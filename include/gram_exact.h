#ifndef GRAM_EXACT_H
#define GRAM_EXACT_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Leonard van Hemert. MIT licensed; see LICENSES/.
 * Exact real Gram matrix lower triangle for row-major binary64 input.
 * A and G must be disjoint, valid allocated array objects; strides in doubles.
 * Every call validates A[i,k]=2^scale_exp X[i,k], integral |X| <=
 * floor_sqrt(floor(2^53/n)); -400<=scale_exp<=200, n<=65536, max_q<=64.
 * For i>0, at most max_q coordinates k<n-1 may differ from A[i-1,k+1].
 * Requires IEEE binary64 correctly rounded operations, round-to-nearest,
 * and compilation without fast-math (use scripts/build.sh).
 * Success 0 guarantees exact lower entries. Upper triangle/padding untouched.
 * m=0 is an unconditional no-op; n=0 writes lower zeros after validation.
 * -1: invalid pointers/sizes/strides; -2: nonfinite/excess defects;
 * -3: allocation failure; -4: unsupported lattice/parameters.
 * On any error discard potentially partial output. No automatic fallback.
 */
int lattice_sequential_gram(const double *A, size_t m, size_t n, size_t lda,
                          double *G, size_t ldg, int scale_exp, size_t max_q);
#ifdef __cplusplus
}
#endif
#endif
