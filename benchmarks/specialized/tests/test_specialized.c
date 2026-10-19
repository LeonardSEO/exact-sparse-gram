#include "hankel_sparse.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
int lattice_sequential_gram(const double*,size_t,size_t,size_t,double*,size_t,int,size_t);
static void compare(const double*A,size_t m,size_t n,size_t lda,int e,size_t q){
 size_t ldg=m+3;double*G=malloc(m*ldg*sizeof(double)),*S=malloc(m*ldg*sizeof(double));assert(G&&S);for(size_t t=0;t<m*ldg;t++)G[t]=S[t]=NAN;
 size_t edits;assert(!hankel_sparse_gram(A,m,n,lda,G,ldg,e,q,&edits));assert(!lattice_sequential_gram(A,m,n,lda,S,ldg,e,q));
 for(size_t i=0;i<m;i++)for(size_t j=0;j<=i;j++){
  int64_t sum=0;for(size_t k=0;k<n;k++)sum+=(int64_t)scalbn(A[i*lda+k],-e)*(int64_t)scalbn(A[j*lda+k],-e);
  assert(G[i*ldg+j]==scalbn((double)sum,2*e));assert(G[i*ldg+j]==S[i*ldg+j]);
 }
 for(size_t i=0;i<m;i++)for(size_t j=i+1;j<ldg;j++)assert(isnan(G[i*ldg+j]));free(G);free(S);
}
int main(void){
 assert(!hankel_sparse_gram(NULL,0,SIZE_MAX,0,NULL,0,999,999,NULL));
 double z[4]={NAN,NAN,NAN,NAN};assert(!hankel_sparse_gram(NULL,2,0,0,z,2,0,0,NULL));assert(z[0]==0&&z[2]==0&&z[3]==0&&isnan(z[1]));
 for(size_t m=1;m<=12;m++)for(size_t n=1;n<=17;n++)for(int e=-400;e<=200;e+=200){
  size_t lda=n+2;double*A=calloc(m*lda,sizeof(double));assert(A);
  for(size_t i=0;i<m;i++)for(size_t k=0;k<n;k++)A[i*lda+k]=scalbn((double)((int)((i*31+k*17)%127)-63),e);
  compare(A,m,n,lda,e,n-1); /* Arbitrary dense E, not just sparse examples. */
  for(size_t i=1;i<m;i++)for(size_t k=0;k+1<n;k++)A[i*lda+k]=A[(i-1)*lda+k+1];
  compare(A,m,n,lda,e,0);free(A);
 }
 for(int e=-400;e<=200;e+=200){
  const int64_t quota=(INT64_C(1)<<53)/7;int64_t b=(int64_t)sqrt((double)quota);while((b+1)*(b+1)<=quota)++b;while(b*b>quota)--b;
  double edge[42];for(size_t i=0;i<6;i++)for(size_t k=0;k<7;k++)edge[i*7+k]=scalbn((double)(((i+k)%2?-1:1)*(b-(int64_t)((i*7+k)%3))),e);compare(edge,6,7,7,e,6);
 }
 double a[6]={0,-67108864.,67108863.,67108864.,-67108864.,0};compare(a,3,2,2,0,1);
 double out[4],bad[4]={0,1,2,3};assert(hankel_sparse_gram(bad,2,2,2,out,2,0,0,NULL)==-2);
 bad[0]=NAN;assert(hankel_sparse_gram(bad,2,2,2,out,2,0,1,NULL)==-2);
 bad[0]=INFINITY;assert(hankel_sparse_gram(bad,2,2,2,out,2,0,1,NULL)==-2);
 bad[0]=0.5;assert(hankel_sparse_gram(bad,2,2,2,out,2,0,1,NULL)==-4);
 bad[0]=nextafter(0.,1.);assert(hankel_sparse_gram(bad,2,2,2,out,2,200,1,NULL)==-4);
 bad[0]=67108865.;assert(hankel_sparse_gram(bad,2,2,2,out,2,0,1,NULL)==-4);
 assert(hankel_sparse_gram(bad,2,2,1,out,2,0,1,NULL)==-1);
 assert(hankel_sparse_gram(bad,2,65537,65537,out,2,0,1,NULL)==-4);
 puts("specialized comparator: full-entry integer oracle, sequential cross-check, padding and admission passed");
}
