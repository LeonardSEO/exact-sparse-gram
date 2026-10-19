/* Fused-prepass SIMD certified adaptive Hankel Gram construction. IEEE binary64,
 * explicit fused fma, correctly rounded sqrt, no fast-math. See notes/numerics.md
 * and notes/communication_guard.md. Independently implemented vector schedule.
 * Lower row-major output; A and G must not overlap. Validation, extraction,
 * row norms, certificates, fresh anchors and temporary allocations are included.
 *
 * Returns 0 success, -1 invalid dimensions/pointers/tolerance, -2 nonfinite or
 * non-Hankel input, -3 allocation failure, -4 outside proved magnitude/n domain,
 * -5 even a direct anchor cannot certify requested tolerance. Except validation
 * errors, an unsuccessful operation may have partially written G.
 */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif

#if FLT_RADIX != 2 || DBL_MANT_DIG != 53
#error "This certificate requires IEEE binary64"
#endif

static double upper_add(double a, double b) {
    if (a == 0.0 && b == 0.0) return 0.0;
    return nextafter(a + b, INFINITY);
}
static double upper_mul(double a, double b) {
    if (a == 0.0 || b == 0.0) return 0.0;
    return nextafter(a * b, INFINITY);
}
static double lower_mul(double a, double b) {
    if (a == 0.0 || b == 0.0) return 0.0;
    return nextafter(a * b, 0.0);
}
static double upper_div(double a, double b) {
    if (a == 0.0) return 0.0;
    return nextafter(a / b, INFINITY);
}
static double gamma_upper(size_t depth) {
    /* depth <= 8195 in the documented domain: integer conversion/scaling exact. */
    const double du = (double)depth * 0x1p-53;
    return upper_div(du, nextafter(1.0 - du, 0.0));
}

#if defined(__aarch64__)
static double reduce8(double q[8]) {
    float64x2_t q0=vld1q_f64(q),q1=vld1q_f64(q+2);
    const float64x2_t q2=vld1q_f64(q+4),q3=vld1q_f64(q+6);
    q0=vaddq_f64(q0,q2); q1=vaddq_f64(q1,q3);
    q0=vaddq_f64(q0,q1);
    return vgetq_lane_f64(q0,0)+vgetq_lane_f64(q0,1);
}

/* All three sums have exactly the same eight chains and binary-tree reduction
 * as scalar dot8; fusion changes memory reuse, not their arithmetic depth. */
static double dot8_fused_norm(const double *a,const double *b,size_t n,
                              double *square,double *absolute) {
    float64x2_t q0=vdupq_n_f64(0),q1=q0,q2=q0,q3=q0;
    float64x2_t s0=q0,s1=q0,s2=q0,s3=q0;
    float64x2_t t0=q0,t1=q0,t2=q0,t3=q0;
    size_t k=0;
    for (; k+8<=n; k+=8) {
#define COMMUNICATION_FUSED_STEP(H) \
        { const float64x2_t av=vld1q_f64(a+k+2*(H)); \
          const float64x2_t bv=vld1q_f64(b+k+2*(H)); \
          q##H=vfmaq_f64(q##H,av,bv); \
          s##H=vfmaq_f64(s##H,vabsq_f64(av),vabsq_f64(bv)); \
          t##H=vfmaq_f64(t##H,av,av); }
        COMMUNICATION_FUSED_STEP(0)
        COMMUNICATION_FUSED_STEP(1)
        COMMUNICATION_FUSED_STEP(2)
        COMMUNICATION_FUSED_STEP(3)
#undef COMMUNICATION_FUSED_STEP
    }
    double q[8],s[8],t[8];
#define COMMUNICATION_STORE(H) \
    vst1q_f64(q+2*(H),q##H); vst1q_f64(s+2*(H),s##H); vst1q_f64(t+2*(H),t##H);
    COMMUNICATION_STORE(0)
    COMMUNICATION_STORE(1)
    COMMUNICATION_STORE(2)
    COMMUNICATION_STORE(3)
#undef COMMUNICATION_STORE
    for (size_t h=0;k<n;++k,++h) {
        q[h]=fma(a[k],b[k],q[h]);
        s[h]=fma(fabs(a[k]),fabs(b[k]),s[h]);
        t[h]=fma(a[k],a[k],t[h]);
    }
    *square=reduce8(t); *absolute=reduce8(s);
    return reduce8(q);
}

static double dot8(const double *a,const double *b,size_t n,double *absolute) {
    float64x2_t q0=vdupq_n_f64(0),q1=q0,q2=q0,q3=q0;
    float64x2_t s0=q0,s1=q0,s2=q0,s3=q0;
    size_t k=0;
    for (; k+8<=n; k+=8) {
#define COMMUNICATION_DOT_STEP(H) \
        { const float64x2_t av=vld1q_f64(a+k+2*(H)); \
          const float64x2_t bv=vld1q_f64(b+k+2*(H)); \
          q##H=vfmaq_f64(q##H,av,bv); \
          s##H=vfmaq_f64(s##H,vabsq_f64(av),vabsq_f64(bv)); }
        COMMUNICATION_DOT_STEP(0)
        COMMUNICATION_DOT_STEP(1)
        COMMUNICATION_DOT_STEP(2)
        COMMUNICATION_DOT_STEP(3)
#undef COMMUNICATION_DOT_STEP
    }
    double q[8],s[8];
#define COMMUNICATION_STORE(H) vst1q_f64(q+2*(H),q##H); vst1q_f64(s+2*(H),s##H);
    COMMUNICATION_STORE(0)
    COMMUNICATION_STORE(1)
    COMMUNICATION_STORE(2)
    COMMUNICATION_STORE(3)
#undef COMMUNICATION_STORE
    for (size_t h=0;k<n;++k,++h) {
        q[h]=fma(a[k],b[k],q[h]);
        s[h]=fma(fabs(a[k]),fabs(b[k]),s[h]);
    }
    if (absolute) *absolute=reduce8(s);
    return reduce8(q);
}

static int validate_input(const double *a,size_t m,size_t n,size_t lda) {
    const float64x2_t zero=vdupq_n_f64(0.0);
    const float64x2_t minimum=vdupq_n_f64(0x1p-400);
    const float64x2_t maximum=vdupq_n_f64(0x1p200);
    for (size_t i=0;i<m;++i) {
        const double *row=a+i*lda;
        const double *previous=i ? row-lda+1 : row;
        const size_t limit=i ? n-1 : n;
        size_t k=0;
        for (;k+8<=limit;k+=8) {
            uint64x2_t valid=vdupq_n_u64(UINT64_MAX);
            for (size_t h=0;h<8;h+=2) {
                const float64x2_t av=vld1q_f64(row+k+h);
                const float64x2_t ab=vabsq_f64(av);
                const uint64x2_t magnitude=vorrq_u64(vceqq_f64(ab,zero),
                    vandq_u64(vcgeq_f64(ab,minimum),vcleq_f64(ab,maximum)));
                valid=vandq_u64(valid,magnitude);
                if (i) valid=vandq_u64(valid,vceqq_f64(av,vld1q_f64(previous+k+h)));
            }
            if (vminvq_u32(vreinterpretq_u32_u64(valid)) != UINT32_MAX) {
                for (size_t h=0;h<8;++h) {
                    const double value=row[k+h];
                    if (!isfinite(value) || (i && value!=previous[k+h])) return -2;
                    if (value!=0.0 && (fabs(value)<0x1p-400 || fabs(value)>0x1p200)) return -4;
                }
            }
        }
        for (;k<n;++k) {
            const double value=row[k];
            if (!isfinite(value) || (i && k+1<n && value!=previous[k])) return -2;
            if (value!=0.0 && (fabs(value)<0x1p-400 || fabs(value)>0x1p200)) return -4;
        }
    }
    return 0;
}
#else
static double dot8(const double *a, const double *b, size_t n, double *absolute) {
    double q[8] = {0}, s[8] = {0};
    size_t k = 0;
    for (; k + 8 <= n; k += 8) {
        for (size_t h = 0; h < 8; ++h) {
            q[h] = fma(a[k+h], b[k+h], q[h]);
            if (absolute) s[h] = fma(fabs(a[k+h]), fabs(b[k+h]), s[h]);
        }
    }
    for (size_t h = 0; k < n; ++k, ++h) {
        q[h] = fma(a[k], b[k], q[h]);
        if (absolute) s[h] = fma(fabs(a[k]), fabs(b[k]), s[h]);
    }
    for (size_t stride = 4; stride; stride /= 2) {
        for (size_t h = 0; h < stride; ++h) {
            q[h] = q[h] + q[h+stride];
            if (absolute) s[h] = s[h] + s[h+stride];
        }
    }
    if (absolute) *absolute = s[0];
    return q[0];
}

static double dot8_fused_norm(const double *a,const double *b,size_t n,
                              double *square,double *absolute) {
    *square=dot8(a,a,n,NULL);
    return dot8(a,b,n,absolute);
}

static int validate_input(const double *A,size_t m,size_t n,size_t lda) {
    for (size_t i=0;i<m;++i) for (size_t k=0;k<n;++k) {
        const double a=A[i*lda+k];
        if (!isfinite(a) || (i && k+1<n && a!=A[(i-1)*lda+k+1])) return -2;
        if (a!=0.0 && (fabs(a)<0x1p-400 || fabs(a)>0x1p200)) return -4;
    }
    return 0;
}
#endif

#if defined(__aarch64__)
/* For nonnegative operands this is nextafter(RN(a*b),0), except that zero
 * remains zero. The representation operation is exact, including infinities
 * and subnormals; no relative-error assumption is used for this helper. */
static float64x2_t lower_mul2(float64x2_t a, float64x2_t b) {
    const float64x2_t zero=vdupq_n_f64(0.0);
    const float64x2_t product=vmulq_f64(a,b);
    const uint64x2_t bits=vreinterpretq_u64_f64(product);
    const float64x2_t predecessor=vreinterpretq_f64_u64(vsubq_u64(bits,vdupq_n_u64(1)));
    return vbslq_f64(vceqq_f64(product,zero),zero,predecessor);
}
#endif

int communication_guard_fast_gram(const double *A, size_t m, size_t n, size_t lda,
                             double *G, size_t ldg, double tolerance,
                             size_t *repairs) {
    if (repairs) *repairs = 0;
    if (m == 0) return 0;
    if (!G || ldg < m || (m-1) > (SIZE_MAX-(m-1))/ldg ||
        !isfinite(tolerance) || tolerance <= 0.0) return -1;
    if (n == 0) {
        for (size_t i=0; i<m; ++i)
            for (size_t j=0; j<=i; ++j) G[i*ldg+j] = 0.0;
        return 0;
    }
    if (!A || lda < n || (m-1) > (SIZE_MAX-(n-1))/lda ||
        m > SIZE_MAX-n+1 || m > SIZE_MAX/sizeof(double)) return -1;
    if (n > 65536) return -4;
    const size_t length=m+n-1;
    if (length > SIZE_MAX/sizeof(double)) return -1;
    const int validation=validate_input(A,m,n,lda);
    if (validation) return validation;
    double *x=malloc(length*sizeof(double));
    double *norm=malloc(m*sizeof(double));
    double *error=malloc(m*sizeof(double));
    if (!x || !norm || !error) {
        free(x); free(norm); free(error); return -3;
    }
    memcpy(x,A,n*sizeof(double));
    for (size_t i=1; i<m; ++i) x[n+i-1]=A[i*lda+n-1];
    const size_t depth=(n+7)/8+3;
    const double gamma=gamma_upper(depth);
    const double anchor_coefficient=upper_div(gamma,nextafter(1.0-gamma,0.0));
    const double norm_denominator=upper_add(1.0,gamma);
    for (size_t i=0; i<m; ++i) {
        double square,absolute;
        G[i*ldg]=dot8_fused_norm(x+i,x,n,&square,&absolute);
        error[i]=upper_mul(anchor_coefficient,absolute);
        if (square==0.0) norm[i]=0.0;
        else {
            const double lower_square=nextafter(square/norm_denominator,0.0);
            norm[i]=nextafter(sqrt(lower_square),0.0);
        }
    }
    int status=0;
    for (size_t i=0; i<m && !status; ++i) {
        double absolute;
        double value=G[i*ldg];
        double bound=error[i];
        double threshold=lower_mul(tolerance,lower_mul(norm[i],norm[0]));
        if (bound>threshold) { status=-5; break; }
        G[i*ldg]=value;
        error[i]=bound;
        size_t j=1;
#if defined(__aarch64__)
        const float64x2_t old_left=vdupq_n_f64(-x[i ? i-1 : 0]);
        const float64x2_t new_left=vdupq_n_f64(x[i+n-1]);
        const float64x2_t eight_u=vdupq_n_f64(0x1p-50);
        const float64x2_t inflate=vdupq_n_f64(0x1.0000000000001p0);
        const float64x2_t tol2=vdupq_n_f64(tolerance);
        const float64x2_t normi=vdupq_n_f64(norm[i]);
        for (; j+1<=i; j+=2) {
            const size_t diagonal=i-j;
            const float64x2_t previous=vld1q_f64(G+(i-1)*ldg+j-1);
            const float64x2_t tmp=vfmaq_f64(previous,old_left,vld1q_f64(x+j-1));
            float64x2_t values=vfmaq_f64(tmp,new_left,vld1q_f64(x+j+n-1));
            float64x2_t previous_error=vld1q_f64(error+diagonal-1);
            previous_error=vextq_f64(previous_error,previous_error,1);
            const float64x2_t sum_abs=vaddq_f64(vabsq_f64(tmp),vabsq_f64(values));
            const float64x2_t increment=vmulq_f64(eight_u,sum_abs);
            float64x2_t bounds=vmulq_f64(vaddq_f64(previous_error,increment),inflate);
            const float64x2_t thresholds=lower_mul2(tol2,lower_mul2(normi,vld1q_f64(norm+j)));
            const float64x2_t maximum=vdupq_n_f64(DBL_MAX);
            const uint64x2_t finite_values=vcleq_f64(vabsq_f64(values),maximum);
            const uint64x2_t finite_bounds=vcleq_f64(bounds,maximum);
            const uint64x2_t failing=vorrq_u64(vcgtq_f64(bounds,thresholds),
                veorq_u64(vandq_u64(finite_values,finite_bounds),vdupq_n_u64(UINT64_MAX)));
            if (vgetq_lane_u64(failing,0) || vgetq_lane_u64(failing,1)) {
                double v[2],e[2],t[2];
                vst1q_f64(v,values); vst1q_f64(e,bounds); vst1q_f64(t,thresholds);
                for (size_t lane=0; lane<2; ++lane) {
                    if (!isfinite(v[lane]) || !isfinite(e[lane]) || e[lane]>t[lane]) {
                        v[lane]=dot8(x+i,x+j+lane,n,&absolute);
                        e[lane]=upper_mul(anchor_coefficient,absolute);
                        if (repairs) ++*repairs;
                        if (e[lane]>t[lane]) { status=-5; break; }
                    }
                }
                if (status) break;
                values=vld1q_f64(v); bounds=vld1q_f64(e);
            }
            vst1q_f64(G+i*ldg+j,values);
            bounds=vextq_f64(bounds,bounds,1);
            vst1q_f64(error+diagonal-1,bounds);
        }
        if (status) break;
#endif
        for (; j<=i; ++j) {
            const size_t diagonal=i-j;
            threshold=lower_mul(tolerance,lower_mul(norm[i],norm[j]));
            const double tmp=fma(-x[i-1],x[j-1],G[(i-1)*ldg+j-1]);
            value=fma(x[i+n-1],x[j+n-1],tmp);
            const double increment=0x1p-50*(fabs(tmp)+fabs(value));
            bound=(error[diagonal]+increment)*0x1.0000000000001p0;
            if (!isfinite(value) || !isfinite(bound) || bound>threshold) {
                value=dot8(x+i,x+j,n,&absolute);
                bound=upper_mul(anchor_coefficient,absolute);
                if (repairs) ++*repairs;
                if (bound>threshold) { status=-5; break; }
            }
            G[i*ldg+j]=value;
            error[diagonal]=bound;
        }
    }
    free(x); free(norm); free(error);
    return status;
}
