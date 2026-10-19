/*
 * Fixed-width real Gram family: feature-major packing and output-register tiles.
 * Untimed communication-track prototype, 2026-10-19.
 *
 * Contract: row-major A[m,n], lower triangle of row-major C[m,m], beta=0.
 * C and workspace must not overlap A or each other. Workspace length is m*n
 * doubles, except n=1 where no workspace is used. The upper triangle is untouched.
 * Caller owns allocation. Time pack+compute together for a fair end-to-end test.
 *
 * Compile: clang -O3 -std=c11 -fno-fast-math -ffp-contract=off -c communication_skinny.c
 * Explicit NEON FMA instructions remain fused independently of fp-contract.
 */
#include <stddef.h>
#include <stdint.h>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif

void communication_pack_features(size_t m, size_t n,
                                 const double *restrict a,
                                 double *restrict packed) {
    for (size_t k = 0; k < n; ++k)
        for (size_t j = 0; j < m; ++j)
            packed[k*m + j] = a[j*n + k];
}

#if defined(__aarch64__)
#define COMMUNICATION_FIXED_KERNEL(N) \
static void communication_kernel_##N(size_t m, const double *restrict a, \
                                      const double *restrict p, double *restrict c) { \
    for (size_t i = 0; i < m; ++i) { \
        const double *ai = a+i*(N); \
        size_t j = 0; \
        for (; j+7 <= i; j += 8) { \
            float64x2_t z0 = vmulq_n_f64(vld1q_f64(p+j), ai[0]); \
            float64x2_t z1 = vmulq_n_f64(vld1q_f64(p+j+2), ai[0]); \
            float64x2_t z2 = vmulq_n_f64(vld1q_f64(p+j+4), ai[0]); \
            float64x2_t z3 = vmulq_n_f64(vld1q_f64(p+j+6), ai[0]); \
            for (size_t k = 1; k < (N); ++k) { \
                const double *pk = p+k*m+j; \
                z0 = vfmaq_n_f64(z0, vld1q_f64(pk), ai[k]); \
                z1 = vfmaq_n_f64(z1, vld1q_f64(pk+2), ai[k]); \
                z2 = vfmaq_n_f64(z2, vld1q_f64(pk+4), ai[k]); \
                z3 = vfmaq_n_f64(z3, vld1q_f64(pk+6), ai[k]); \
            } \
            vst1q_f64(c+i*m+j, z0); \
            vst1q_f64(c+i*m+j+2, z1); \
            vst1q_f64(c+i*m+j+4, z2); \
            vst1q_f64(c+i*m+j+6, z3); \
        } \
        for (; j <= i; ++j) { \
            double z = ai[0]*p[j]; \
            for (size_t k = 1; k < (N); ++k) z = __builtin_fma(ai[k], p[k*m+j], z); \
            c[i*m+j] = z; \
        } \
    } \
}
#else
#define COMMUNICATION_FIXED_KERNEL(N) \
static void communication_kernel_##N(size_t m, const double *restrict a, \
                                      const double *restrict p, double *restrict c) { \
    for (size_t i = 0; i < m; ++i) \
        for (size_t j = 0; j <= i; ++j) { \
            double z = a[i*(N)]*p[j]; \
            for (size_t k = 1; k < (N); ++k) \
                z = __builtin_fma(a[i*(N)+k], p[k*m+j], z); \
            c[i*m+j] = z; \
        } \
}
#endif

COMMUNICATION_FIXED_KERNEL(1)
COMMUNICATION_FIXED_KERNEL(2)
COMMUNICATION_FIXED_KERNEL(3)
COMMUNICATION_FIXED_KERNEL(4)
COMMUNICATION_FIXED_KERNEL(6)
COMMUNICATION_FIXED_KERNEL(8)
COMMUNICATION_FIXED_KERNEL(12)
COMMUNICATION_FIXED_KERNEL(16)
COMMUNICATION_FIXED_KERNEL(24)
COMMUNICATION_FIXED_KERNEL(32)

/* Returns zero on supported dimensions; -1 for unsupported feature counts. */
int communication_gram_prepared(size_t m, size_t n, const double *restrict a,
                                const double *restrict p, double *restrict c) {
    switch (n) {
    case 0:
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j <= i; ++j) c[i*m+j] = 0.0;
        return 0;
#define COMMUNICATION_CASE(N) case N: communication_kernel_##N(m,a,p,c); return 0;
    COMMUNICATION_CASE(1)
    COMMUNICATION_CASE(2)
    COMMUNICATION_CASE(3)
    COMMUNICATION_CASE(4)
    COMMUNICATION_CASE(6)
    COMMUNICATION_CASE(8)
    COMMUNICATION_CASE(12)
    COMMUNICATION_CASE(16)
    COMMUNICATION_CASE(24)
    COMMUNICATION_CASE(32)
#undef COMMUNICATION_CASE
    default: return -1;
    }
}

int communication_gram(size_t m, size_t n, const double *restrict a,
                       double *restrict workspace, double *restrict c) {
    if (n == 1) return communication_gram_prepared(m, n, a, a, c);
    switch (n) {
    case 0: case 2: case 3: case 4: case 6: case 8: case 12:
    case 16: case 24: case 32: break;
    default: return -1;
    }
    communication_pack_features(m,n,a,workspace);
    return communication_gram_prepared(m,n,a,workspace,c);
}
