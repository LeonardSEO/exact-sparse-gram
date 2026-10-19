#ifndef HANKEL_SPARSE_H
#define HANKEL_SPARSE_H
#include <stddef.h>
/* Same lower-output domain/status/stride contract as lattice_sequential_gram.
 * stats, if nonnull, receives actual nnz(A-H) on success. */
int hankel_sparse_gram(const double*,size_t,size_t,size_t,double*,size_t,int,size_t,size_t*);
#endif
