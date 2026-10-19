#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fenv.h>
#include <errno.h>
#include <limits.h>
#include <inttypes.h>
#include "hankel_sparse.h"
int lattice_sequential_gram(const double*,size_t,size_t,size_t,double*,size_t,int,size_t);
static uint64_t state;
static uint64_t rnd(void){state^=state<<13;state^=state>>7;state^=state<<17;return state;}
static volatile double consume;
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC_RAW,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void lattice_mirror(double *G,size_t m){
 const size_t b=32;
 for(size_t bi=0;bi<m;bi+=b)for(size_t bj=bi;bj<m;bj+=b){
  size_t ie=bi+b<m?bi+b:m,je=bj+b<m?bj+b:m;
  for(size_t i=bi;i<ie;i++)for(size_t j=bj>i?bj:i+1;j<je;j++)G[i*m+j]=G[j*m+i];
 }
}
static int lattice_call(int method,const double*A,size_t m,size_t n,double*G,int e,size_t q,int full){
 int status=0;
 if(method==0)cblas_dsyrk(CblasRowMajor,CblasLower,CblasNoTrans,(int)m,(int)n,1.,A,(int)n,0.,G,(int)m);
 else if(method==4)cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)m,(int)m,(int)n,1.,A,(int)n,A,(int)n,0.,G,(int)m);
 else if(method==13)status=lattice_sequential_gram(A,m,n,n,G,m,e,q);
 else if(method==14)status=hankel_sparse_gram(A,m,n,n,G,m,e,q,NULL);
 else return -10;
 if(!status&&full&&method!=4)lattice_mirror(G,m);return status;
}
static double coefficient(int dist,size_t k){
 if(dist==0)return (double)((int)(rnd()&65535)-32768);
 if(dist==1)return 16384.+(double)((int)(rnd()%17)-8);
 if(dist==2)return (k&1)?-32768.:32767.;
 if(dist==3)return k%17?0.:((k&1)?-32768.:32767.);
 return (double)((int)(rnd()%3)-1);
}
static void lattice_make(double*A,size_t m,size_t n,int accumulating,size_t q,int e,int dist,uint64_t seed){
 state=seed?seed:1;double scale=scalbn(1.,e);
 if(accumulating){
  for(size_t k=0;k<n;k++)A[k]=coefficient(dist,k)*scale;
  for(size_t i=1;i<m;i++){
   memcpy(A+i*n,A+(i-1)*n+1,(n-1)*sizeof(double));
   A[i*n+n-1]=coefficient(dist,i+n-1)*scale;
   for(size_t h=0;h<q;h++){
    size_t k=(size_t)(rnd()%(n-1));
    /* Repeated positions only lower actualq. Each modification is validated. */
    A[i*n+k]=coefficient(dist,i*n+k)*scale;
   }
  }
 }else{
  double*x=malloc((m+n-1)*sizeof(double));
  for(size_t k=0;k<m+n-1;k++)x[k]=coefficient(dist,k)*scale;
  for(size_t i=0;i<m;i++)memcpy(A+i*n,x+i,n*sizeof(double));
  free(x);
  if(q){
   /* One independent edit per window guarantees q<=2. */
   for(size_t i=0;i<m;i++)A[i*n+(17*i+7)%n]=0.;
  }
 }
}
static int64_t icoeff(double x,int e){return (int64_t)scalbn(x,-e);}
static int verify_samples(const double*A,const double*G,size_t m,size_t n,int e,int method,uint64_t seed,int full){
 /* Independent scalar int64products/sums. Coefficient dot bounds <=2^46. */
 size_t ndiag=m<128?m:128;size_t count=ndiag+256;
 uint64_t saved=state;state=seed+991;
 for(size_t s=0;s<count;s++){
  size_t i=s<ndiag?s:(size_t)(rnd()%m),j=s<ndiag?s:(size_t)(rnd()%m);
  if(!full&&i<j){size_t v=i;i=j;j=v;}
  int64_t sum=0;
  for(size_t k=0;k<n;k++)sum+=icoeff(A[i*n+k],e)*icoeff(A[j*n+k],e);
  double exact=scalbn((double)sum,2*e);
  double got=G[i*m+j];
  if(method==13||method==14){if(got!=exact){state=saved;return -1;}}
  else {
   /* Vendor/oldcertificate commonpresetnormscaled accuracy, not exactness. */
   double ni=0,nj=0;for(size_t k=0;k<n;k++){ni=fma(A[i*n+k],A[i*n+k],ni);nj=fma(A[j*n+k],A[j*n+k],nj);}
   if(!isfinite(got)||fabs(got-exact)>2e-12*sqrt(ni)*sqrt(nj)){state=saved;return -1;}
  }
 }
 state=saved;return 0;
}
static int parse_unsigned(const char*s,uint64_t*out){
 char*end;errno=0;if(!*s||*s=='-'||*s=='+')return -1;unsigned long long v=strtoull(s,&end,10);if(errno||*end)return -1;*out=(uint64_t)v;return 0;
}
static int parse_int(const char*s,int*out){
 char*end;errno=0;if(!*s)return -1;long v=strtol(s,&end,10);if(errno||*end||v<INT_MIN||v>INT_MAX)return -1;*out=(int)v;return 0;
}
/* Canonical little-endian bytes of binary64 bit patterns, no RNG consumption. */
static uint64_t fingerprint(const double*A,size_t count){uint64_t hash=UINT64_C(14695981039346656037);for(size_t i=0;i<count;i++){uint64_t bits;memcpy(&bits,A+i,8);for(size_t b=0;b<8;b++){hash^=(bits>>(8*b))&255;hash*=UINT64_C(1099511628211);}}return hash;}
static void make_orders(uint64_t seed,int samples,int orders[1001][4]){
 uint64_t saved=state;state=seed+713;int base[4];const int methods[4]={0,4,13,14};
 for(int s=0;s<samples;s++){
  if(s%4==0){for(int z=0;z<4;z++)base[z]=z;for(int z=3;z>0;z--){int j=(int)(rnd()%(unsigned)(z+1)),v=base[z];base[z]=base[j];base[j]=v;}}
  for(int z=0;z<4;z++)orders[s][z]=methods[base[(z+s%4)%4]];
 }state=saved;
}
static void print_orders(int samples,int orders[1001][4]){printf("[ ");for(int s=0;s<samples;s++)printf("%s[%d,%d,%d,%d]",s?",":"",orders[s][0],orders[s][1],orders[s][2],orders[s][3]);printf("]");}
int main(int argc,char**argv){
 int diagnostic=argc==12&&!strcmp(argv[11],"--input-diagnostic");
 if(argc!=11&&!diagnostic){fprintf(stderr,"usage: specialized_benchmark windows|innov m n q e seed dist samples full batch_ms [--input-diagnostic]\n");return 2;}
 uint64_t um,un,uq,seed;int e,dist,samples,full;
 if(parse_unsigned(argv[2],&um)||parse_unsigned(argv[3],&un)||parse_unsigned(argv[4],&uq)||parse_int(argv[5],&e)||parse_unsigned(argv[6],&seed)||parse_int(argv[7],&dist)||parse_int(argv[8],&samples)||parse_int(argv[9],&full))return 2;
 if(um>SIZE_MAX||un>SIZE_MAX||uq>SIZE_MAX)return 2;
 size_t m=(size_t)um,n=(size_t)un,q=(size_t)uq;
 char*end;errno=0;double target=strtod(argv[10],&end)*.001;if(errno||!*argv[10]||*end)return 2;
 if(strcmp(argv[1],"innov")&&strcmp(argv[1],"windows"))return 2;
 if(!m||!n||m>4096||n>65536||q>64||samples<4||samples>1001||samples%4||e<-200||e>180||dist<0||dist>4||(full!=0&&full!=1)||!isfinite(target)||target<=0)return 2;
 if((!strcmp(argv[1],"innov")&&n==1&&q)||(!strcmp(argv[1],"windows")&&q!=0&&q!=2))return 2;
 if(fegetround()!=FE_TONEAREST||BLASSetThreading(BLAS_THREADING_SINGLE_THREADED)||BLASGetThreading()!=BLAS_THREADING_SINGLE_THREADED)return 3;
 double*A=malloc(m*n*sizeof(double)),*G=malloc(m*m*sizeof(double));if(!A||!G)return 4;
 lattice_make(A,m,n,!strcmp(argv[1],"innov"),q,e,dist,seed);
 size_t actual_q=0,max_observed_q=0;
 for(size_t i=1;i<m;i++){
  size_t row_q=0;for(size_t k=0;k+1<n;k++)row_q+=A[i*n+k]!=A[(i-1)*n+k+1];
  actual_q+=row_q;if(row_q>max_observed_q)max_observed_q=row_q;
 }
 uint64_t input_hash=fingerprint(A,m*n);int orders[1001][4];make_orders(seed,samples,orders);
 if(diagnostic){printf("{\"family\":\"%s\",\"m\":%zu,\"n\":%zu,\"q_cap\":%zu,\"scale_exp\":%d,\"seed\":%" PRIu64 ",\"dist\":%d,\"output_full\":%d,\"input_fingerprint_fnv1a64\":\"%016" PRIx64 "\",\"generator_final_state\":%" PRIu64 ",\"method_orders\":",argv[1],m,n,q,e,seed,dist,full,input_hash,state);print_orders(samples,orders);puts("}");free(A);free(G);return 0;}
 const int methods[4]={0,4,13,14};
 size_t actual_e=0; if(hankel_sparse_gram(A,m,n,n,G,m,e,q,&actual_e))return 5;
 for(int z=0;z<4;z++){
  /* Untimed poison prevents another method's earlier output hiding a miss. */
  for(size_t t=0;t<m*m;t++)G[t]=NAN;
  for(int w=0;w<5;w++){
   int status=lattice_call(methods[z],A,m,n,G,e,q,full);
   if(status){fprintf(stderr,"method%d status%d\n",methods[z],status);return 5;}
  }
  for(size_t i=0;i<m;i++)for(size_t j=0;j<(full?m:i+1);j++)if(!isfinite(G[i*m+j]))return 6;
  if(verify_samples(A,G,m,n,e,methods[z],seed,full)){fprintf(stderr,"exactoracle failure method%d\n",methods[z]);return 6;}
  if(full&&methods[z]!=4)for(size_t i=0;i<m;i++)for(size_t j=0;j<i;j++)if(G[i*m+j]!=G[j*m+i])return 6;
 }
 size_t batch[4];
 for(int z=0;z<4;z++){
  size_t b=1;double elapsed;
  do {double t=now();for(size_t k=0;k<b;k++)if(lattice_call(methods[z],A,m,n,G,e,q,full))return 7;elapsed=now()-t;if(elapsed<target)b*=2;}while(elapsed<target);
  batch[z]=b;
 }
 printf("{\"family\":\"%s\",\"m\":%zu,\"n\":%zu,\"q_cap\":%zu,\"scale_exp\":%d,\"seed\":%llu,\"dist\":%d,\"samples\":[",argv[1],m,n,q,e,(unsigned long long)seed,dist);
 for(int s=0;s<samples;s++){
  double times[4];
  for(int zz=0;zz<4;zz++){
   int z=0;while(methods[z]!=orders[s][zz])++z;double t=now();
   for(size_t k=0;k<batch[z];k++)if(lattice_call(methods[z],A,m,n,G,e,q,full))return 7;
   times[z]=(now()-t)/batch[z];consume+=G[(m-1)*m];
  }
  printf("%s{",s?",":"");for(int z=0;z<4;z++)printf("%s\"%d\":%.17g",z?",":"",methods[z],times[z]);printf("}");fflush(stdout);
 }
 printf("],\"input_fingerprint_fnv1a64\":\"%016" PRIx64 "\",\"actual_e\":%zu,\"output_full\":%d,\"actual_q\":%zu,\"max_observed_q\":%zu,\"blas_threads\":%d,\"oracle\":\"independentint64samplechecks+allfiniteafterpoison\",\"batches\":{",input_hash,actual_e,full,actual_q,max_observed_q,BLASGetThreading());
 for(int z=0;z<4;z++)printf("%s\"%d\":%zu",z?",":"",methods[z],batch[z]);printf("},\"method_orders\":");print_orders(samples,orders);puts("}");
 free(A);free(G);return 0;
}
