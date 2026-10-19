/* Certified adaptive Hankel Gram construction. IEEE binary64, round-to-nearest,
 * explicit fused fma, correctly rounded sqrt, no fast-math. See notes/numerics.md.
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

static int guard_impl(const double *A, size_t m, size_t n, size_t lda,
                      double *G, size_t ldg, size_t anchor, double tolerance,
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
        for (size_t j=0; j<=i; ++j) {
            const size_t diagonal=i-j;
            const double threshold=lower_mul(tolerance,lower_mul(norm[i],norm[j]));
            double value, bound;
            const int forced=(j==0 || (anchor && j%anchor==0));
            if (!forced) {
                const double tmp=fma(-x[i-1],x[j-1],G[(i-1)*ldg+j-1]);
                value=fma(x[i+n-1],x[j+n-1],tmp);
                /* 4u*fl(|tmp|+|value|) >= (u/(1-u))*(|tmp|+|value|).
                 * Multiplication by 4u is an exact power-of-two scaling here. */
                const double increment=0x1p-51*(fabs(tmp)+fabs(value));
                bound=upper_add(error[diagonal],increment);
            } else { value=0.0; bound=INFINITY; }
            if (forced || bound>threshold) {
                double absolute;
                value=dot8(x+i,x+j,n,&absolute);
                bound=upper_mul(anchor_coefficient,absolute);
                if (!forced && repairs) ++*repairs;
                if (bound>threshold) { status=-5; break; }
            }
            G[i*ldg+j]=value;
            error[diagonal]=bound;
        }
    }
    free(x); free(norm); free(error);
    return status;
}

int numerical_guard_gram(const double *A,size_t m,size_t n,size_t lda,
                         double *G,size_t ldg,size_t anchor,double tolerance) {
    return guard_impl(A,m,n,lda,G,ldg,anchor,tolerance,NULL);
}

int algebra_hankel_gram_guard(const double *A,size_t m,size_t n,size_t lda,
                             double *G,size_t ldg,double tolerance,size_t *repairs) {
    return guard_impl(A,m,n,lda,G,ldg,0,tolerance,repairs);
}
