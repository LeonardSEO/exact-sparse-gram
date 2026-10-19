/* Exact binary64 Gram by sequential factor updates on sparse row shifts.
 * Derived from this project's src/lattice_substitute.c, SHA-256:
 * d1400c003886f0e57102e75bb7e442b7329cbad26d6d796b844b9bf877f1ba95.
 * The validation/anchor scaffold is project-local; no third-party source or
 * new license grant is introduced. See notes/lattice_sequential.md.
 * API: row-major A/G, strides in doubles, disjoint valid array objects;
 * write lower G only. Every call validates A[i,k] = 2^e X[i,k], X integral,
 * |X| <= floor_sqrt(floor(2^53/n)); e in [-400,200], 1 <= n <= 65536.
 * Each shifted row has <= max_q defects, max_q <= 64. No quantization.
 * Sequentially replace the left factors, then replace the right factors;
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
#if defined(__aarch64__) && !defined(LATTICE_SEQUENTIAL_SCALAR)
#include <arm_neon.h>
#define LATTICE_NEON 1
#endif

_Static_assert(DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 &&
               DBL_MIN_EXP == -1021 && sizeof(double) == 8,
               "lattice_sequential_gram requires IEEE binary64 double");

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

/* Only used after a block predicate fails; preserve scalar error classes.
 * Reaching the final return would mean inconsistent scalar/vector predicates;
 * reject conservatively rather than accepting without completing the scan.
 */
static int failed_block_status(const double *a, double scale,
                               double inverse, double bound) {
    for (size_t h = 0; h < 8; ++h) {
        const int status = scalar_valid(a[h], scale, inverse, bound);
        if (status) return status;
    }
    return -4;
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
    /* Strict inequality reserves a scalar tail of 1..8 elements. Every
     * shifted vector ends at previous[k+8], which is inside the row. */
    for (; k + 8 < n; k += 8) {
        float64x2_t av[4];
        uint64x2_t bounded = vdupq_n_u64(UINT64_MAX);
#if defined(__clang__)
#pragma clang loop unroll(full)
#endif
        for (size_t h = 0; h < 4; ++h) {
            av[h] = vld1q_f64(a + k + 2 * h);
            bounded = vandq_u64(bounded, vcleq_f64(vabsq_f64(av[h]), bv));
        }
        /* The finite bound rejects NaN and infinities too. One branch
         * validates the range of all eight values. */
        if (!both_lanes(bounded))
            return failed_block_status(a + k, scale, inverse, bound);
        uint64x2_t on_grid = vdupq_n_u64(UINT64_MAX);
#if defined(__clang__)
#pragma clang loop unroll(full)
#endif
        for (size_t h = 0; h < 4; ++h) {
            const float64x2_t x = vmulq_f64(av[h], iv);
            const uint64x2_t integral = vceqq_f64(x, vrndq_f64(x));
            const uint64x2_t reconstructed = vceqq_f64(vmulq_f64(x, sv), av[h]);
            on_grid = vandq_u64(on_grid, vandq_u64(integral, reconstructed));
        }
        /* One branch checks integrality AND reconstruction for all eight. */
        if (!both_lanes(on_grid))
            return failed_block_status(a + k, scale, inverse, bound);
        if (previous) {
            uint64x2_t same = vdupq_n_u64(UINT64_MAX);
#if defined(__clang__)
#pragma clang loop unroll(full)
#endif
            for (size_t h = 0; h < 4; ++h)
                same = vandq_u64(same, vceqq_f64(av[h],
                                    vld1q_f64(previous + k + 2 * h + 1)));
            /* One sparse branch per eight features. Extraction stays in
             * increasing coordinate order, including every changed lane. */
            if (!both_lanes(same)) {
                for (size_t h = 0; h < 8; ++h) {
                    const size_t c = k + h;
                    const int status = append_residual(a[c], previous[c + 1], c,
                                                        capq, &count, used, r);
                    if (status) return status;
                }
            }
        }
#if defined(__clang__)
#pragma clang loop unroll(full)
#endif
        for (size_t h = 0; h < 4; ++h)
            qv[h] = vfmaq_f64(qv[h], av[h], vld1q_f64(first + k + 2 * h));
    }
    for (size_t h = 0; h < 4; ++h) vst1q_f64(q + 2 * h, qv[h]);
#else
    for (; k + 8 < n; k += 8) {
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

int lattice_sequential_gram(const double *A, size_t m, size_t n, size_t lda,
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
                    /* First group replaces Ui[k] by Ai[k], with every
                     * right factor still Uj[k]. Second replaces Uj[k] by
                     * Aj[k], with every left factor already Ai[k]. Each
                     * prefix is a sum of n bounded mixed-factor products. */
                    for (size_t t = rows[i]; t < rows[i + 1]; ++t)
                        value = fma(r[t].d, uj[r[t].k + 1], value);
                    for (size_t t = rows[j]; t < rows[j + 1]; ++t)
                        value = fma(ai[r[t].k], r[t].d, value);
                }
                G[i * ldg + j] = value;
            }
        }
    }
    free(r); free(rows);
    return status;
}
