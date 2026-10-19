#include <stddef.h>
#include <stdlib.h>
int communication_gram(size_t,size_t,const double*,double*,double*);
int communication_alloc_gram(const double*A,size_t m,size_t n,double*G){
 double *P=n==1?NULL:malloc(m*n*sizeof(double));
 if(n!=1&&!P)return -2;
 int ret=communication_gram(m,n,A,P,G);free(P);return ret;
}
