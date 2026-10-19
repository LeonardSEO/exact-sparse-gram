/* SIMD certified adaptive Hankel Gram construction. IEEE binary64, round-to-nearest,
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

int communication_guard_gram(const double *A, size_t m, size_t n, size_t lda,
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
    for (size_t i=0; i<m; ++i) {
        for (size_t k=0; k<n; ++k) {
            double a=A[i*lda+k];
            if (!isfinite(a)) return -2;
            if (i && k+1<n && a != A[(i-1)*lda+k+1]) return -2;
            if (a != 0.0 && (fabs(a)<0x1p-400 || fabs(a)>0x1p200)) return -4;
        }
    }
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
        const double square=dot8(x+i,x+i,n,NULL);
        if (square==0.0) norm[i]=0.0;
        else {
            const double lower_square=nextafter(square/norm_denominator,0.0);
            norm[i]=nextafter(sqrt(lower_square),0.0);
        }
    }
    int status=0;
    for (size_t i=0; i<m && !status; ++i) {
        double absolute;
        double value=dot8(x+i,x,n,&absolute);
        double bound=upper_mul(anchor_coefficient,absolute);
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
