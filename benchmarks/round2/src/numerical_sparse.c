/* Adaptive certified sparse row-shift Gram algorithm. binary64 round-to-nearest,
 * explicit fma, correctly rounded sqrt; no fast-math; disjoint A/G arrays.
 * Sparse residuals use exact TwoDiff hi+lo, never rounded residuals alone.
 * Domain: n<=65536; max_q<=64; each nonzero |A_ik| in [2^-400,2^200].
 * Same error target as numerical_guard.c: |error_ij|<=tol*norm(row_i)*norm(row_j).
 * Codes 0 success; -1 sizes/tolerance; -2 nonfinite or too many residuals;
 * -3 allocation; -4 outside proved magnitude/n/q domain; -5 certificate fails.
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

struct numerical_residual { size_t k; double hi,lo; };
static double upadd(double a,double b) {
    return (a==0.0 && b==0.0)?0.0:nextafter(a+b,INFINITY);
}
static double upmul(double a,double b) {
    return (a==0.0 || b==0.0)?0.0:nextafter(a*b,INFINITY);
}
static double downmul(double a,double b) {
    return (a==0.0 || b==0.0)?0.0:nextafter(a*b,0.0);
}
static double updiv(double a,double b) {
    return a==0.0?0.0:nextafter(a/b,INFINITY);
}
static void two_diff(double a,double b,double *hi,double *lo) {
    const double x=a-b, bv=a-x, av=x+bv, br=bv-b, ar=a-av;
    *hi=x; *lo=ar+br;
}
static double sparse_dot8(const double *a,const double *b,size_t n,double *absolute) {
    double q[8]={0},s[8]={0}; size_t k=0;
    for (; k+8<=n; k+=8) {
        for (size_t h=0; h<8; ++h) {
            q[h]=fma(a[k+h],b[k+h],q[h]);
            if (absolute) s[h]=fma(fabs(a[k+h]),fabs(b[k+h]),s[h]);
        }
    }
    for (size_t h=0; k<n; ++k,++h) {
        q[h]=fma(a[k],b[k],q[h]);
        if (absolute) s[h]=fma(fabs(a[k]),fabs(b[k]),s[h]);
    }
    for (size_t stride=4; stride; stride/=2) {
        for (size_t h=0; h<stride; ++h) {
            q[h]+=q[h+stride];
            if (absolute) s[h]+=s[h+stride];
        }
    }
    if (absolute) *absolute=s[0];
    return q[0];
}

int numerical_sparse_guard_gram(const double *A,size_t m,size_t n,size_t lda,
                                double *G,size_t ldg,size_t anchor,size_t max_q,
                                double tolerance,size_t *residual_count,
                                size_t *repairs) {
    if (residual_count) *residual_count=0;
    if (repairs) *repairs=0;
    if (m==0) return 0;
    if (!G || ldg<m || (m-1)>(SIZE_MAX-(m-1))/ldg ||
        !isfinite(tolerance) || tolerance<=0.0) return -1;
    if (n==0) {
        for (size_t i=0;i<m;++i)
            for (size_t j=0;j<=i;++j) G[i*ldg+j]=0.0;
        return 0;
    }
    if (!A || lda<n || (m-1)>(SIZE_MAX-(n-1))/lda || m==SIZE_MAX)
        return -1;
    if (n>65536 || max_q>64) return -4;
    const size_t capq=max_q<n-1?max_q:n-1;
    if (capq && (m-1)>SIZE_MAX/capq) return -1;
    const size_t capacity=(m-1)*capq;
    if (capacity>SIZE_MAX/sizeof(struct numerical_residual) ||
        m+1>SIZE_MAX/sizeof(size_t) || m>SIZE_MAX/sizeof(double)) return -1;
    size_t *rows=malloc((m+1)*sizeof(size_t));
    struct numerical_residual *r=capacity?malloc(capacity*sizeof(*r)):NULL;
    double *norm=malloc(m*sizeof(double)), *error=malloc(m*sizeof(double));
    if (!rows || (capacity && !r) || !norm || !error) {
        free(rows);free(r);free(norm);free(error);return -3;
    }
    int status=0; size_t used=0;
    for (size_t i=0;i<m && !status;++i) {
        rows[i]=used; size_t count=0;
        for (size_t k=0;k<n;++k) {
            const double a=A[i*lda+k];
            if (!isfinite(a)) {status=-2;break;}
            if (a!=0.0 && (fabs(a)<0x1p-400 || fabs(a)>0x1p200)) {status=-4;break;}
            if (i && k+1<n && a!=A[(i-1)*lda+k+1]) {
                if (count==capq) {status=-2;break;}
                struct numerical_residual *v=&r[used++];v->k=k;
                two_diff(a,A[(i-1)*lda+k+1],&v->hi,&v->lo);
                ++count;
            }
        }
        rows[i+1]=used;
    }
    if (status) {free(rows);free(r);free(norm);free(error);return status;}
    if (residual_count) *residual_count=used;
    const size_t depth=(n+7)/8+3;
    const double du=(double)depth*0x1p-53;
    const double gamma=updiv(du,nextafter(1.0-du,0.0));
    const double coefficient=updiv(gamma,nextafter(1.0-gamma,0.0));
    const double denominator=upadd(1.0,gamma);
    for (size_t i=0;i<m;++i) {
        const double square=sparse_dot8(A+i*lda,A+i*lda,n,NULL);
        norm[i]=square==0.0?0.0:
            nextafter(sqrt(nextafter(square/denominator,0.0)),0.0);
    }
    for (size_t i=0;i<m && !status;++i) {
        for (size_t j=0;j<=i;++j) {
            const size_t diagonal=i-j;
            const double threshold=downmul(tolerance,downmul(norm[i],norm[j]));
            const int forced=(j==0 || (anchor && j%anchor==0));
            double value=0.0,bound=INFINITY;
            if (!forced) {
                value=fma(-A[(i-1)*lda],A[(j-1)*lda],G[(i-1)*ldg+j-1]);
                double magnitude=fabs(value);
                value=fma(A[i*lda+n-1],A[j*lda+n-1],value);
                magnitude+=fabs(value);
                for (size_t t=rows[i];t<rows[i+1];++t) {
                    const double a=A[j*lda+r[t].k];
                    value=fma(r[t].hi,a,value);magnitude+=fabs(value);
                    if (r[t].lo!=0.0) {
                        value=fma(r[t].lo,a,value);magnitude+=fabs(value);
                    }
                }
                for (size_t t=rows[j];t<rows[j+1];++t) {
                    const double a=A[(i-1)*lda+r[t].k+1];
                    value=fma(a,r[t].hi,value);magnitude+=fabs(value);
                    if (r[t].lo!=0.0) {
                        value=fma(a,r[t].lo,value);magnitude+=fabs(value);
                    }
                }
                /* <=258 FMA outputs: 4u*fl(sum|out|) safely exceeds
                 * rho*sum|out|; one outward addition preserves old bound. */
                bound=upadd(error[diagonal],0x1p-51*magnitude);
            }
            if (forced || bound>threshold) {
                double absolute;
                value=sparse_dot8(A+i*lda,A+j*lda,n,&absolute);
                bound=upmul(coefficient,absolute);
                if (!forced && repairs) ++*repairs;
                if (bound>threshold) {status=-5;break;}
            }
            G[i*ldg+j]=value;error[diagonal]=bound;
        }
    }
    free(rows);free(r);free(norm);free(error);return status;
}

int algebra_sparse_gram_guard(const double *A,size_t m,size_t n,size_t lda,
                              double *G,size_t ldg,double tolerance,size_t max_q,
                              size_t *residual_count,size_t *repairs) {
    return numerical_sparse_guard_gram(A,m,n,lda,G,ldg,0,max_q,tolerance,
                                      residual_count,repairs);
}
