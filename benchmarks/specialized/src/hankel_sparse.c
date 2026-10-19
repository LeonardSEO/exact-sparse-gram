/* Independent project implementation of established sliding dots + H+E.
 * No third-party implementation was copied. See docs/SPECIALIZED_COMPARATOR.md.
 * Strict IEEE binary64 arithmetic; no fast-math. */
#include "hankel_sparse.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif
_Static_assert(DBL_MANT_DIG==53 && DBL_MAX_EXP==1024 && DBL_MIN_EXP==-1021 && sizeof(double)==8,"requires binary64");
struct edit { size_t k; double d; };
static int extent(size_t r,size_t c,size_t s){size_t lim=SIZE_MAX/sizeof(double);return !r||!c||(c<=lim&&s>=c&&r-1<=(lim-c)/s);}
static uint64_t bound_for(size_t n){uint64_t q=(UINT64_C(1)<<53)/n,b=(uint64_t)sqrt((double)q);while((b+1)*(b+1)<=q)++b;while(b*b>q)--b;return b;}
static int valid(double a,double scale,double inv,double bound){if(!isfinite(a))return -2;if(fabs(a)>bound)return -4;double x=a*inv;return x==trunc(x)&&x*scale==a?0:-4;}
/* Batch range/grid predicates, with scalar fallback retaining error classes. */
static int validate_row(const double*a,size_t n,double scale,double inv,double bound){
 size_t k=0;
#if defined(__aarch64__)
 float64x2_t sv=vdupq_n_f64(scale),iv=vdupq_n_f64(inv),bv=vdupq_n_f64(bound);
 for(;k+8<=n;k+=8){uint64x2_t ok=vdupq_n_u64(UINT64_MAX);
  for(size_t h=0;h<4;h++){float64x2_t v=vld1q_f64(a+k+2*h),x=vmulq_f64(v,iv);ok=vandq_u64(ok,vandq_u64(vcleq_f64(vabsq_f64(v),bv),vandq_u64(vceqq_f64(x,vrndq_f64(x)),vceqq_f64(vmulq_f64(x,sv),v))));}
  if(!(vgetq_lane_u64(ok,0)&vgetq_lane_u64(ok,1))){for(size_t h=0;h<8;h++){int s=valid(a[k+h],scale,inv,bound);if(s)return s;}return -4;}
 }
#endif
 for(;k<n;k++){int s=valid(a[k],scale,inv,bound);if(s)return s;}return 0;
}
static double dot(const double*a,const double*b,size_t n){double v[8]={0};size_t k=0;for(;k+8<=n;k+=8)for(size_t h=0;h<8;h++)v[h]=fma(a[k+h],b[k+h],v[h]);for(size_t h=0;k<n;k++,h++)v[h]=fma(a[k],b[k],v[h]);for(size_t s=4;s;s/=2)for(size_t h=0;h<s;h++)v[h]+=v[h+s];return v[0];}
int hankel_sparse_gram(const double*A,size_t m,size_t n,size_t lda,double*G,size_t ldg,int e,size_t maxq,size_t*stats){
 if(stats)*stats=0;
 if(!m)return 0;
 if(!G||!extent(m,m,ldg))return -1;
 if(e< -400||e>200||n>65536||maxq>64)return -4;
 if(!n){for(size_t i=0;i<m;i++)for(size_t j=0;j<=i;j++)G[i*ldg+j]=0;return 0;}
 if(!A||!extent(m,n,lda)||m==SIZE_MAX||m>SIZE_MAX-n+1)return -1;
 size_t hn=m+n-1;if(hn>SIZE_MAX/sizeof(double)||(m+1)>SIZE_MAX/sizeof(size_t))return -1;
 double*h=malloc(hn*sizeof(*h));size_t*rows=malloc((m+1)*sizeof(*rows));struct edit*r=NULL;
 if(!h||!rows){free(h);free(rows);return -3;}
 for(size_t k=0;k<n;k++)h[k]=A[k];
 for(size_t i=1;i<m;i++)h[n+i-1]=A[i*lda+n-1];
 size_t used=0;int status=0;double scale=scalbn(1,e),inv=scalbn(1,-e),bound=(double)bound_for(n)*scale;
 /* Detect actual E and local transition defects from dense A; nothing latent. */
 for(size_t i=0;i<m;i++){
  const double*a=A+i*lda;rows[i]=used;status=validate_row(a,n,scale,inv,bound);if(status)goto done;
  size_t q=0;for(size_t k=0;k<n;k++){
   if(a[k]!=h[i+k]){if(used==SIZE_MAX){status=-1;goto done;}++used;}
   if(i&&k+1<n&&a[k]!=A[(i-1)*lda+k+1]){if(++q>maxq){status=-2;goto done;}}
  }
 }
 rows[m]=used;
 if(used>SIZE_MAX/sizeof(*r)){status=-1;goto done;}
 if(used){r=malloc(used*sizeof(*r));if(!r){status=-3;goto done;}}
 for(size_t i=0;i<m;i++){size_t t=rows[i];for(size_t k=0;k<n;k++)if(A[i*lda+k]!=h[i+k]){r[t].k=k;r[t++].d=A[i*lda+k]-h[i+k];}}
 /* H-only phase must finish before in-place corrections overwrite HH^T. */
 for(size_t i=0;i<m;i++){
  G[i*ldg]=dot(h+i,h,n);
  for(size_t j=1;j<=i;j++){double v=fma(-h[i-1],h[j-1],G[(i-1)*ldg+j-1]);G[i*ldg+j]=fma(h[i+n-1],h[j+n-1],v);}
 }
 for(size_t i=0;i<m;i++)for(size_t j=0;j<=i;j++){
  double v=G[i*ldg+j];
  for(size_t t=rows[i];t<rows[i+1];t++)v=fma(r[t].d,h[j+r[t].k],v);
  for(size_t t=rows[j];t<rows[j+1];t++)v=fma(A[i*lda+r[t].k],r[t].d,v);
  G[i*ldg+j]=v;
 }
 if(stats)*stats=used;
 done:free(r);free(rows);free(h);return status;
}
