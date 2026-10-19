#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/resource.h>
#include <fenv.h>

int joint_gram_lower(const double*,size_t,size_t,size_t,double*,size_t);
int algebra_hankel_gram(const double*,size_t,size_t,size_t,double*,size_t,size_t);
int algebra_hankel_gram_guard(const double*,size_t,size_t,size_t,double*,size_t,double,size_t*);
int communication_alloc_gram(const double*,size_t,size_t,double*);
int communication_guard_gram(const double*,size_t,size_t,size_t,double*,size_t,double,size_t*);
int algebra_sparse_gram(const double*,size_t,size_t,size_t,double*,size_t,size_t,size_t,size_t*);
int algebra_sparse_gram_guard(const double*,size_t,size_t,size_t,double*,size_t,double,size_t,size_t*,size_t*);
int communication_guard_fast_gram(const double*,size_t,size_t,size_t,double*,size_t,double,size_t*);
int numerical_sparse_fast_gram(const double*,size_t,size_t,size_t,double*,size_t,double,size_t,size_t*,size_t*);
static int output_full;
static volatile double consume;
static uint64_t state;
static uint64_t rnd(void){state^=state<<13;state^=state>>7;state^=state<<17;return state;}
static double uni(void){return (double)(rnd()>>11)*0x1p-53*2.-1.;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC_RAW,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void syrk(const double*A,size_t m,size_t n,double*G){
 cblas_dsyrk(CblasRowMajor,CblasLower,CblasNoTrans,(int)m,(int)n,1.,A,(int)n,0.,G,(int)m);
}
static void gemm(const double*A,size_t m,size_t n,double*G){
 cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)m,(int)m,(int)n,1.,A,(int)n,A,(int)n,0.,G,(int)m);
}
static int call_core(int method,const double*A,size_t m,size_t n,double*G,size_t anchor){
 if(method==0){syrk(A,m,n,G);return 0;}
 if(method==1)return joint_gram_lower(A,m,n,n,G,m);
 if(method==5)return communication_alloc_gram(A,m,n,G);
 if(method==2)return algebra_hankel_gram(A,m,n,n,G,m,anchor);
#ifndef SKIP_GUARD
 if(method==3)return algebra_hankel_gram_guard(A,m,n,n,G,m,2e-12,NULL);
#endif
 if(method==9)return numerical_sparse_fast_gram(A,m,n,n,G,m,2e-12,2,NULL,NULL);
 if(method==10)return communication_guard_fast_gram(A,m,n,n,G,m,2e-12,NULL);
 if(method==6)return communication_guard_gram(A,m,n,n,G,m,2e-12,NULL);
 if(method==7)return algebra_sparse_gram(A,m,n,n,G,m,anchor,2,NULL);
 if(method==8)return algebra_sparse_gram_guard(A,m,n,n,G,m,2e-12,2,NULL,NULL);
 gemm(A,m,n,G);return 0;
}
static int call(int method,const double*A,size_t m,size_t n,double*G,size_t anchor){
 int status=call_core(method,A,m,n,G,anchor);
 if(!status && output_full && method!=4){
  const size_t block=32;
  for(size_t bi=0;bi<m;bi+=block)for(size_t bj=bi;bj<m;bj+=block){
   size_t ie=bi+block<m?bi+block:m,je=bj+block<m?bj+block:m;
   for(size_t i=bi;i<ie;i++)for(size_t j=(bj>i?bj:i+1);j<je;j++)G[i*m+j]=G[j*m+i];
  }
 }
 return status;
}
static void make(double*A,size_t m,size_t n,int hankel,int dist,uint64_t seed){
 state=seed?seed:1;
 if(hankel){
  double*x=malloc((m+n-1)*sizeof(double));
  for(size_t k=0;k<m+n-1;k++){
   if(dist==0)x[k]=uni();
   else if(dist==1){x[k]=0;for(int q=0;q<12;q++)x[k]+=uni()/2;}
   else if(dist==2)x[k]=1.+0x1p-24*uni();
   else if(dist==3)x[k]=(k%17==0?1.:0.)*(k%2?-1.:1.);
   else if(dist==4)x[k]=ldexp(uni(),(int)(rnd()%401)-200);
   else x[k]=(k%2?-1.:1.)*(1.+0x1p-30*uni());
  }
  for(size_t i=0;i<m;i++)for(size_t k=0;k<n;k++)A[i*n+k]=x[i+k];
  free(x);
 }else for(size_t k=0;k<m*n;k++){
  if(dist==0)A[k]=uni();
  else if(dist==1){A[k]=0;for(int q=0;q<12;q++)A[k]+=uni()/2;}
  else if(dist==2)A[k]=1.+0x1p-24*uni();
  else if(dist==3)A[k]=(k%17==0?1.:0.)*(k%2?-1.:1.);
  else if(dist==4)A[k]=ldexp(uni(),(int)(rnd()%401)-200);
  else A[k]=(k%2?-1.:1.)*(1.+0x1p-30*uni());
 }
}
int main(int argc,char**argv){
 if(argc<7){fprintf(stderr,"usage: benchmark dense|hankel m n seed dist samples [anchor] [batch_ms]\n");return 2;}
 size_t m=strtoull(argv[2],0,10),n=strtoull(argv[3],0,10),anchor=argc>7?strtoull(argv[7],0,10):32;
 uint64_t seed=strtoull(argv[4],0,10);int hankel=strcmp(argv[1],"dense")!=0,dist=atoi(argv[5]),samples=atoi(argv[6]);
 double target=argc>8?atof(argv[8])*.001:.010;
 if(!m||!n||m>16384||n>1048576||samples<1||samples>1001)return 2;
 if(fegetround()!=FE_TONEAREST){fprintf(stderr,"requires round-to-nearest\n");return 3;}
 int set_result=BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);
 if(set_result || BLASGetThreading()!=BLAS_THREADING_SINGLE_THREADED){fprintf(stderr,"BLAS single-thread setting failed %d %d\n",set_result,BLASGetThreading());return 3;}
 double*A=malloc(m*n*8),*G=malloc(m*m*8),*R=malloc(m*m*8);
 if(!A||!G||!R)return 4;
 make(A,m,n,hankel,dist,seed);
 if(!strcmp(argv[1],"sparse")){for(size_t i=0;i<m;i++)A[i*n+((i*17+7)%n)]=0.;}
 output_full=getenv("GRAM_FULL") && atoi(getenv("GRAM_FULL"));
 syrk(A,m,n,R);
 int methods[6]={0,4,1,5},count=4;
 if(hankel){methods[2]=2;methods[3]=10;count=4;
#ifdef SKIP_GUARD
 count=3;
#endif
 }
 if(!strcmp(argv[1],"sparse")){methods[2]=8;methods[3]=9;count=4;}
 for(int z=0;z<count;z++){
  int method=methods[z];for(int w=0;w<5;w++)if(call(method,A,m,n,G,anchor)!=0){fprintf(stderr,"method %d failed\n",method);return 5;}
  double worst=0.;
  for(size_t i=0;i<m;i++)for(size_t j=0;j<(output_full?m:i+1);j++){
   double denom=sqrt(R[i*m+i])*sqrt(R[j*m+j]);double reference=i>=j?R[i*m+j]:R[j*m+i];double err=fabs(G[i*m+j]-reference);
   if(!isfinite(G[i*m+j]) || (denom==0&&err!=0)){fprintf(stderr,"nonfinite/zero-row failure\n");return 6;}
   if(denom>0&&err/denom>worst)worst=err/denom;
  }
  fprintf(stderr,"check method=%d m=%zu n=%zu relative_norm_vs_syrk=%.6g\n",method,m,n,worst);
  if(worst>4e-12)return 6;
 }
 size_t batches[6]={0};
 for(int z=0;z<count;z++){
  int method=methods[z];size_t batch=1;double elapsed;
  do {double t=now();for(size_t b=0;b<batch;b++)call(method,A,m,n,G,anchor);elapsed=now()-t;if(elapsed<target)batch*=2;}while(elapsed<target);
  batches[z]=batch;
 }
 printf("{\"mode\":\"%s\",\"m\":%zu,\"n\":%zu,\"seed\":%llu,\"dist\":%d,\"anchor\":%zu,\"blas_threads\":%d,\"samples\":[",argv[1],m,n,(unsigned long long)seed,dist,anchor,BLASGetThreading());
 state=seed+713;
 int base_order[6];
 for(int s=0;s<samples;s++){
  int order[6];
  if(s%count==0){for(int z=0;z<count;z++)base_order[z]=z;for(int z=count-1;z>0;z--){int q=(int)(rnd()%(z+1)),v=base_order[z];base_order[z]=base_order[q];base_order[q]=v;}}
  for(int z=0;z<count;z++)order[z]=base_order[(z+s%count)%count];
  double times[6];
  for(int zz=0;zz<count;zz++){
   int z=order[zz],method=methods[z];double t=now();
   for(size_t b=0;b<batches[z];b++)if(call(method,A,m,n,G,anchor)!=0)return 7;
   times[z]=(now()-t)/batches[z];consume+=G[(m-1)*m];
  }
  if(s)printf(",");printf("{");
  for(int z=0;z<count;z++)printf("%s\"%d\":%.17g",z?",":"",methods[z],times[z]);
  printf("}");fflush(stdout);
 }
 printf("],\"output_full\":%d,\"batches\":{",output_full);for(int z=0;z<count;z++)printf("%s\"%d\":%zu",z?",":"",methods[z],batches[z]);
 printf("}}\n");free(A);free(G);free(R);return 0;
}
