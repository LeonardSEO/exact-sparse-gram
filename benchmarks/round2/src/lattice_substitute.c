/* Exact binary64 Gram by coordinate substitution on sparse row shifts.
 * Derived from this project's src/lattice_exact.c, SHA-256:
 * 3402d7c6531fdb17178493df4205b54aaa580fcc6f6d99aa3a66cffa3e3da540.
 * The validation/anchor scaffold is project-local; no third-party source or
 * new license grant is introduced. See notes/lattice_substitute.md.
 * API: row-major A/G, strides in doubles, disjoint valid array objects;
 * write lower G only. Every call validates A[i,k] = 2^e X[i,k], X integral,
 * |X| <= floor_sqrt(floor(2^53/n)); e in [-400,200], 1 <= n <= 65536.
 * Each shifted row has <= max_q defects, max_q <= 64. No quantization.
 * Merging defect coordinates replaces old coordinate products by new ones;
 * every accumulator has <= n bounded products, independent of q and m.
 * Compile without fast-math. IEEE binary64 correctly rounded +,-,*,fma required.
 * Return 0 success, -1 invalid sizes/pointers, -2 nonfinite/too many defects,
 * -3 allocation failure, -4 outside the supported exact lattice domain.
 * m=0 is a no-op. n=0 writes lower zeros after size/parameter validation.
 * On error the caller must discard potentially partial lower output.
 */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#if defined(__aarch64__) && !defined(LATTICE_SUBSTITUTE_SCALAR)
#include <arm_neon.h>
#define LATTICE_NEON 1
#endif

_Static_assert(DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 &&
               DBL_MIN_EXP == -1021 && sizeof(double) == 8,
               "lattice_substitute_gram requires IEEE binary64 double");

struct lattice_residual { size_t k; double d; };

static int valid_extent(size_t rows, size_t cols, size_t stride) {
    const size_t limit = SIZE_MAX / sizeof(double);
    if (!rows || !cols) return 1;
    return cols <= limit && stride >= cols &&
           rows - 1 <= (limit - cols) / stride;
}

/* n is nonzero and <= 65536. sqrt supplies an initial estimate only;
 * uint64 arithmetic proves the final floor bound. All squared candidates are
 * below 2^54, far below uint64 overflow. No input coefficients are trusted.
 */
static uint64_t coefficient_bound(size_t n) {
    const uint64_t quota = (UINT64_C(1) << 53) / (uint64_t)n;
    uint64_t b = (uint64_t)sqrt((double)quota);
    while ((b + 1) * (b + 1) <= quota) ++b;
    while (b * b > quota) --b;
    return b;
}

static int scalar_valid(double a, double scale, double inverse, double bound) {
    if (!isfinite(a)) return -2;
    if (fabs(a) > bound) return -4;
    const double x = a * inverse;
    /* Reconstruction is necessary: an invalid tiny a can scale to zero. */
    return x == trunc(x) && x * scale == a ? 0 : -4;
}

static int append_residual(double a, double b, size_t k, size_t capq,
                           size_t *count, size_t *used,
                           struct lattice_residual *r) {
    if (a == b) return 0;
    if (*count == capq) return -2;
    r[*used].k = k;
    r[*used].d = a - b; /* Exact: integer coefficient magnitude <= 2B < 2^28. */
    ++*used;
    ++*count;
    return 0;
}

#if defined(LATTICE_NEON)
static uint64_t both_lanes(uint64x2_t mask) {
    return vgetq_lane_u64(mask, 0) & vgetq_lane_u64(mask, 1);
}

static int vector_valid(float64x2_t a, float64x2_t scale,
                        float64x2_t inverse, float64x2_t bound) {
    const float64x2_t ab = vabsq_f64(a);
    if (!both_lanes(vcleq_f64(ab, vdupq_n_f64(DBL_MAX)))) return -2;
    if (!both_lanes(vcleq_f64(ab, bound))) return -4;
    const float64x2_t x = vmulq_f64(a, inverse);
    const uint64x2_t integral = vceqq_f64(x, vrndq_f64(x));
    const uint64x2_t reconstructed = vceqq_f64(vmulq_f64(x, scale), a);
    return both_lanes(vandq_u64(integral, reconstructed)) ? 0 : -4;
}
#endif

/* One feature traversal validates the lattice, extracts sparse residuals, and
 * computes a first-column anchor using eight independent exact FMA chains.
 * The predecessor row and row zero were validated before their use here.
 */
static int validate_seed(const double *a, const double *first,
                         const double *previous, size_t n, size_t capq,
                         double scale, double inverse, double bound,
                         struct lattice_residual *r, size_t *used,
                         double *anchor) {
    double q[8] = {0};
    size_t k = 0, count = 0;
#if defined(LATTICE_NEON)
    const float64x2_t sv = vdupq_n_f64(scale), iv = vdupq_n_f64(inverse);
    const float64x2_t bv = vdupq_n_f64(bound);
    float64x2_t qv[4] = {vdupq_n_f64(0), vdupq_n_f64(0),
                         vdupq_n_f64(0), vdupq_n_f64(0)};
    for (; k + 8 <= n; k += 8) {
#if defined(__clang__)
#pragma clang loop unroll(full)
#endif
        for (size_t h = 0; h < 4; ++h) {
            const size_t c = k + 2 * h;
            const float64x2_t av = vld1q_f64(a + c);
            int status = vector_valid(av, sv, iv, bv);
            if (status) return status;
            if (previous) {
                if (c + 2 < n) {
                    const float64x2_t pv = vld1q_f64(previous + c + 1);
                    if (!both_lanes(vceqq_f64(av, pv))) {
                        status = append_residual(a[c], previous[c + 1], c,
                                                  capq, &count, used, r);
                        if (status) return status;
                        status = append_residual(a[c + 1], previous[c + 2], c + 1,
                                                  capq, &count, used, r);
                        if (status) return status;
                    }
                } else {
                    /* This vector ends at column n-1: only its first lane
                     * has a shifted predecessor. Never read previous[n]. */
                    status = append_residual(a[c], previous[c + 1], c,
                                              capq, &count, used, r);
                    if (status) return status;
                }
            }
            qv[h] = vfmaq_f64(qv[h], av, vld1q_f64(first + c));
        }
    }
    for (size_t h = 0; h < 4; ++h) vst1q_f64(q + 2 * h, qv[h]);
#else
    for (; k + 8 <= n; k += 8) {
        for (size_t h = 0; h < 8; ++h) {
            const size_t c = k + h;
            int status = scalar_valid(a[c], scale, inverse, bound);
            if (status) return status;
            if (previous && c + 1 < n) {
                status = append_residual(a[c], previous[c + 1], c,
                                          capq, &count, used, r);
                if (status) return status;
            }
            q[h] = fma(a[c], first[c], q[h]);
        }
    }
#endif
    for (size_t h = 0; k < n; ++k, ++h) {
        int status = scalar_valid(a[k], scale, inverse, bound);
        if (status) return status;
        if (previous && k + 1 < n) {
            status = append_residual(a[k], previous[k + 1], k,
                                      capq, &count, used, r);
            if (status) return status;
        }
        q[h] = fma(a[k], first[k], q[h]);
    }
    for (size_t stride = 4; stride; stride /= 2)
        for (size_t h = 0; h < stride; ++h) q[h] += q[h + stride];
    *anchor = q[0];
    return 0;
}

int lattice_substitute_gram(const double *A, size_t m, size_t n, size_t lda,
                       double *G, size_t ldg, int scale_exp, size_t max_q) {
    if (!m) return 0;
    if (!G || !valid_extent(m, m, ldg)) return -1;
    if (scale_exp < -400 || scale_exp > 200 || n > 65536 || max_q > 64)
        return -4;
    if (!n) {
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j <= i; ++j) G[i * ldg + j] = 0.0;
        return 0;
    }
    if (!A || !valid_extent(m, n, lda) || m == SIZE_MAX) return -1;
    const size_t capq = max_q < n - 1 ? max_q : n - 1;
    if (capq && m - 1 > SIZE_MAX / capq) return -1;
    const size_t capacity = (m - 1) * capq;
    if (capacity > SIZE_MAX / sizeof(struct lattice_residual) ||
        (capacity && m + 1 > SIZE_MAX / sizeof(size_t))) return -1;
    struct lattice_residual *r = capacity ? malloc(capacity * sizeof(*r)) : NULL;
    size_t *rows = capacity ? malloc((m + 1) * sizeof(*rows)) : NULL;
    if (capacity && (!r || !rows)) {
        free(r); free(rows); return -3;
    }
    const double scale = scalbn(1.0, scale_exp);
    const double inverse = scalbn(1.0, -scale_exp);
    const uint64_t b = coefficient_bound(n);
    const double bound = (double)b * scale;
    size_t used = 0;
    int status = 0;
    for (size_t i = 0; i < m; ++i) {
        if (rows) rows[i] = used;
        status = validate_seed(A + i * lda, A,
                                i ? A + (i - 1) * lda : NULL,
                                n, capq, scale, inverse, bound, r, &used,
                                G + i * ldg);
        if (status) break;
    }
    if (!status && rows) rows[m] = used;
    if (!status) {
        for (size_t i = 1; i < m; ++i) {
            const double *ai = A + i * lda;
            const double *ui = A + (i - 1) * lda;
            for (size_t j = 1; j <= i; ++j) {
                const double *aj = A + j * lda;
                const double *uj = A + (j - 1) * lda;
                double value = fma(-ui[0], uj[0],
                                    G[(i - 1) * ldg + j - 1]);
                value = fma(ai[n - 1], aj[n - 1], value);
                if (rows) {
                    size_t ti = rows[i], tj = rows[j];
                    const size_t ei = rows[i + 1], ej = rows[j + 1];
                    while (ti < ei || tj < ej) {
                        if (tj == ej || (ti < ei && r[ti].k < r[tj].k)) {
                            /* Dj[k]=0: old Ui*Uj -> new Ai*Aj. */
                            value = fma(r[ti].d, aj[r[ti].k], value);
                            ++ti;
                        } else if (ti == ei || r[tj].k < r[ti].k) {
                            /* Di[k]=0: old Ui*Uj -> new Ai*Aj. */
                            value = fma(ui[r[tj].k + 1], r[tj].d, value);
                            ++tj;
                        } else {
                            /* Both changed: delete old, then insert new.
                             * There are n-1 bounded products after deletion,
                             * and n after insertion; never the grouped cross
                             * term responsible for the extra 2q B^2 bound. */
                            const size_t k = r[ti].k;
                            value = fma(-ui[k + 1], uj[k + 1], value);
                            value = fma(ai[k], aj[k], value);
                            ++ti;
                            ++tj;
                        }
                    }
                }
                G[i * ldg + j] = value;
            }
        }
    }
    free(r); free(rows);
    return status;
}
