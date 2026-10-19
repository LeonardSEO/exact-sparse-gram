/* Sparse row-shift displacement Gram construction. Research artifact.
 * Exact identity: let U_i[k]=A[i-1,k+1], D_i[k]=A[i,k]-U_i[k] (k<n-1).
 * G_ij=G_(i-1,j-1)-A[i-1,0]A[j-1,0]+A[i,n-1]A[j,n-1]
 *      +D_i dot A_j[0:n-1]+U_i dot D_j.
 * D is represented as exact TwoDiff hi+lo, assuming IEEE binary64,
 * round-to-nearest, gradual underflow, and no intermediate overflow.
 * cc -O3 -std=c11 -fno-fast-math -c algebra_sparse.c
 * Lower triangle only; output/preallocated input row-major and may have padding.
 * A and G must refer to nonoverlapping storage, as for ordinary BLAS SYRK.
 * anchor0: seed first column only; b>=1: fresh seed when j%b==0.
 * max_q bounds the number of nonzero D_i entries per input row.
 * Return0 success,-1 invalid dimensions,-2 nonfinite/out-of-family/overflow,
 * -3 allocation failure. G untouched on validation/allocation failure.
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

struct algebra_residual { size_t k; double hi, lo; };

static void algebra_two_diff(double a, double b, double *hi, double *lo) {
    const double x = a - b;
    const double bv = a - x;
    const double av = x + bv;
    const double br = bv - b;
    const double ar = a - av;
    *hi = x;
    *lo = ar + br;
}

int algebra_sparse_gram(const double *A, size_t m, size_t n, size_t lda,
                        double *G, size_t ldg, size_t anchor, size_t max_q,
                        size_t *residual_count) {
    if (residual_count) *residual_count = 0;
    if (m == 0) return 0;
    if (!G || ldg < m || (m - 1) > (SIZE_MAX - (m - 1)) / ldg) return -1;
    if (n == 0) {
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j <= i; ++j) G[i * ldg + j] = 0.0;
        return 0;
    }
    if (!A || lda < n || (m - 1) > (SIZE_MAX - (n - 1)) / lda || m == SIZE_MAX)
        return -1;
    const size_t capq = max_q < n - 1 ? max_q : n - 1;
    if (capq != 0 && (m - 1) > SIZE_MAX / capq) return -1;
    const size_t capacity = (m - 1) * capq;
    if (capacity > SIZE_MAX / sizeof(struct algebra_residual) ||
        m + 1 > SIZE_MAX / sizeof(size_t)) return -1;
    size_t *rows = (size_t *)malloc((m + 1) * sizeof(size_t));
    struct algebra_residual *r = capacity ?
        (struct algebra_residual *)malloc(capacity * sizeof(*r)) : NULL;
    if (!rows || (capacity && !r)) { free(rows); free(r); return -3; }
    size_t used = 0;
    rows[0] = 0;
    for (size_t i = 0; i < m; ++i) {
        rows[i] = used;
        size_t rowq = 0;
        for (size_t k = 0; k < n; ++k) {
            const double a = A[i * lda + k];
            if (!isfinite(a)) { free(rows); free(r); return -2; }
            if (i != 0 && k + 1 < n && a != A[(i - 1) * lda + k + 1]) {
                if (rowq == capq) { free(rows); free(r); return -2; }
                struct algebra_residual *v = &r[used++];
                v->k = k;
                algebra_two_diff(a, A[(i - 1) * lda + k + 1], &v->hi, &v->lo);
                if (!isfinite(v->hi) || !isfinite(v->lo)) {
                    free(rows); free(r); return -2;
                }
                ++rowq;
            }
        }
        rows[i + 1] = used;
    }
    if (residual_count) *residual_count = used;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double q;
            if (j == 0 || (anchor != 0 && j % anchor == 0)) {
                q = 0.0;
                for (size_t k = 0; k < n; ++k)
                    q = fma(A[i * lda + k], A[j * lda + k], q);
            } else {
                q = fma(-A[(i - 1) * lda], A[(j - 1) * lda],
                        G[(i - 1) * ldg + j - 1]);
                q = fma(A[i * lda + n - 1], A[j * lda + n - 1], q);
                for (size_t t = rows[i]; t < rows[i + 1]; ++t) {
                    const double a = A[j * lda + r[t].k];
                    q = fma(r[t].hi, a, q);
                    if (r[t].lo != 0.0) q = fma(r[t].lo, a, q);
                }
                for (size_t t = rows[j]; t < rows[j + 1]; ++t) {
                    const double a = A[(i - 1) * lda + r[t].k + 1];
                    q = fma(a, r[t].hi, q);
                    if (r[t].lo != 0.0) q = fma(a, r[t].lo, q);
                }
            }
            if (!isfinite(q)) { free(rows); free(r); return -2; }
            G[i * ldg + j] = q;
        }
    }
    free(rows);
    free(r);
    return 0;
}
