/* Exact-real Hankel Gram construction with periodic fresh anchors.
 * Copyright: research artifact, 2026-10-19. No claim of novelty for recurrence.
 * Build: cc -O3 -std=c11 -fno-fast-math -c algebra_hankel.c
 * ABI: row-major A/G; computes only G[i*ldg+j], 0 <= j <= i < m.
 * anchor == 0: one initial seed per diagonal; anchor >= 1: fresh dot at
 * each column j divisible by anchor. anchor == 1 computes every dot fresh.
 * Returns 0 success; -1 invalid sizes/pointers; -2 nonfinite/non-Hankel A;
 * -3 allocation failure. On -2 G is untouched.
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int algebra_hankel_gram(const double *A, size_t m, size_t n, size_t lda,
                        double *G, size_t ldg, size_t anchor) {
    if (m == 0) return 0;
    if (!G || ldg < m || (m - 1) > (SIZE_MAX - (m - 1)) / ldg) return -1;
    if (n == 0) {
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j <= i; ++j) G[i * ldg + j] = 0.0;
        return 0;
    }
    if (!A || lda < n || (m - 1) > (SIZE_MAX - (n - 1)) / lda || m > SIZE_MAX - n + 1)
        return -1;
    const size_t length = m + (n - 1);
    if (length > SIZE_MAX / sizeof(double)) return -1;
    /* Validation is part of every operation, not trusted preprocessing. */
    for (size_t i = 0; i < m; ++i) {
        for (size_t k = 0; k < n; ++k) {
            const double a = A[i * lda + k];
            if (!isfinite(a)) return -2;
            if (i != 0 && k + 1 < n && a != A[(i - 1) * lda + k + 1])
                return -2;
        }
    }
    double *x = (double *)malloc(length * sizeof(double));
    if (!x) return -3;
    memcpy(x, A, n * sizeof(double));
    for (size_t i = 1; i < m; ++i) x[n + i - 1] = A[i * lda + n - 1];
    /* Row-major traversal minimizes output traffic. Previous row is live. */
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double q;
            if (j == 0 || (anchor != 0 && j % anchor == 0)) {
                q = 0.0;
                for (size_t k = 0; k < n; ++k)
                    q = fma(x[i + k], x[j + k], q);
            } else {
                q = fma(-x[i - 1], x[j - 1], G[(i - 1) * ldg + j - 1]);
                q = fma(x[i + n - 1], x[j + n - 1], q);
            }
            G[i * ldg + j] = q;
        }
    }
    free(x);
    return 0;
}
