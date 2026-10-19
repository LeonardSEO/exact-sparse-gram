/*
 * Feature-major packed, register-tiled Gram kernels.
 *
 * Contract: A and G are non-overlapping, row-major binary64 arrays;
 * lda >= n and ldg >= m.  Write only G[i*ldg+j], 0 <= j <= i < m.
 * Return 0 on success, -1 on allocation/size failure, -2 on invalid strides.
 * All input-dependent packing/allocation/free is inside the public call.
 * Compile: clang -O3 -std=c11 -fPIC -shared src/joint_kernel.c -o libjoint.dylib
 * Do not enable -ffast-math when using the stated error bound.
 */
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#if defined(__aarch64__) || defined(__arm64__)
#include <arm_neon.h>
#define JOINT_NEON 1
#else
#define JOINT_NEON 0
#endif

static int checked_product(size_t a, size_t b, size_t *out) {
    if (b && a > SIZE_MAX / b) return -1;
    *out = a * b;
    return 0;
}

static double *pack_features(const double *restrict A, size_t m, size_t n,
                             size_t lda, size_t *stride) {
    size_t count, bytes;
    if (m > SIZE_MAX - 7) return NULL;
    *stride = (m + 7) & ~(size_t)7;
    if (checked_product(n, *stride, &count) ||
        checked_product(count, sizeof(double), &bytes)) return NULL;
    double *P = (double *)malloc(bytes);
    if (!P) return NULL;
    /* k outside keeps writes contiguous; for skinny n A remains cache-resident. */
    for (size_t k = 0; k < n; ++k)
        for (size_t i = 0; i < m; ++i)
            P[k * *stride + i] = A[i * lda + k];
    return P;
}

#if JOINT_NEON

/* Four output rows, eight consecutive columns, entirely below the diagonal.
 * 16 vector accumulators are independent; each feature loads 12 input doubles
 * and performs 32 binary64 fused multiply-adds (64 useful real operations).
 */
static inline void tile4x8(const double *restrict P, size_t ps, size_t n,
                          size_t i, size_t j, double *restrict G, size_t gs) {
    const double *p = P;
    float64x2_t x0 = vld1q_f64(p+j), x1 = vld1q_f64(p+j+2);
    float64x2_t x2 = vld1q_f64(p+j+4), x3 = vld1q_f64(p+j+6);
#define INIT(R) \
    float64x2_t c##R##0 = vmulq_n_f64(x0,p[i+R]); \
    float64x2_t c##R##1 = vmulq_n_f64(x1,p[i+R]); \
    float64x2_t c##R##2 = vmulq_n_f64(x2,p[i+R]); \
    float64x2_t c##R##3 = vmulq_n_f64(x3,p[i+R]);
    INIT(0) INIT(1) INIT(2) INIT(3)
#undef INIT
    for (size_t k = 1; k < n; ++k) {
        p = P + k * ps;
        x0 = vld1q_f64(p+j);   x1 = vld1q_f64(p+j+2);
        x2 = vld1q_f64(p+j+4); x3 = vld1q_f64(p+j+6);
#define UPDATE(R) \
        c##R##0 = vfmaq_n_f64(c##R##0,x0,p[i+R]); \
        c##R##1 = vfmaq_n_f64(c##R##1,x1,p[i+R]); \
        c##R##2 = vfmaq_n_f64(c##R##2,x2,p[i+R]); \
        c##R##3 = vfmaq_n_f64(c##R##3,x3,p[i+R]);
        UPDATE(0) UPDATE(1) UPDATE(2) UPDATE(3)
#undef UPDATE
    }
#define STORE(R) \
    vst1q_f64(G+(i+R)*gs+j,c##R##0); \
    vst1q_f64(G+(i+R)*gs+j+2,c##R##1); \
    vst1q_f64(G+(i+R)*gs+j+4,c##R##2); \
    vst1q_f64(G+(i+R)*gs+j+6,c##R##3);
    STORE(0) STORE(1) STORE(2) STORE(3)
#undef STORE
}

static inline void tail_row(const double *restrict P, size_t ps, size_t n,
                            size_t i, size_t j, double *restrict G, size_t gs) {
    for (; j + 2 <= i + 1; j += 2) {
        float64x2_t c = vmulq_n_f64(vld1q_f64(P+j),P[i]);
        for (size_t k = 1; k < n; ++k) {
            const double *p = P + k * ps;
            c = vfmaq_n_f64(c,vld1q_f64(p+j),p[i]);
        }
        vst1q_f64(G+i*gs+j,c);
    }
    if (j <= i) {
        double c = P[i] * P[j];
        for (size_t k = 1; k < n; ++k)
            c = __builtin_fma(P[k*ps+i],P[k*ps+j],c);
        G[i*gs+j] = c;
    }
}

#else
static inline void tail_row(const double *restrict P, size_t ps, size_t n,
                            size_t i, size_t j, double *restrict G, size_t gs) {
    for (; j <= i; ++j) {
        double c = P[i] * P[j];
        for (size_t k = 1; k < n; ++k)
            c = __builtin_fma(P[k*ps+i],P[k*ps+j],c);
        G[i*gs+j] = c;
    }
}
#endif

/* Workspace worker: P[k*ps+i] = A[i*lda+k], ps >= m; n > 0.
 * Public primarily to allow independent exact-verification harnesses.
 */
int joint_gram_lower_packed(const double *P, size_t m, size_t n, size_t ps,
                           double *G, size_t ldg) {
    if (ps < m || ldg < m || !n) return -2;
    size_t i = 0;
#if JOINT_NEON
    for (; i + 4 <= m; i += 4) {
        size_t j = 0;
        for (; j + 8 <= i; j += 8)
            tile4x8(P,ps,n,i,j,G,ldg);
        tail_row(P,ps,n,i,  j,G,ldg);
        tail_row(P,ps,n,i+1,j,G,ldg);
        tail_row(P,ps,n,i+2,j,G,ldg);
        tail_row(P,ps,n,i+3,j,G,ldg);
    }
#endif
    for (; i < m; ++i) tail_row(P,ps,n,i,0,G,ldg);
    return 0;
}

int joint_gram_lower(const double *A, size_t m, size_t n, size_t lda,
                     double *G, size_t ldg) {
    if (lda < n || ldg < m) return -2;
    if (!m) return 0;
    if (!n) {
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j <= i; ++j) G[i*ldg+j] = 0.0;
        return 0;
    }
    /* n=1 needs no conversion: rows already form one contiguous feature. */
    if (n == 1 && lda == 1)
        return joint_gram_lower_packed(A,m,n,m,G,ldg);
    size_t ps;
    double *P = pack_features(A,m,n,lda,&ps);
    if (!P) return -1;
    int status = joint_gram_lower_packed(P,m,n,ps,G,ldg);
    free(P);
    return status;
}
