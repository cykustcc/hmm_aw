#include <math.h>

#define __USE_C99_MATH
#include "utils/blas_like.h"
#include "utils/blas_utils.h"
#include "math.h"
#include <stdio.h>
#include <assert.h>

void _dgzero(size_t n, double *a) {
  size_t i;
  for (i=0; i<n; ++i) assert(a[i] > 1E-10);
}

void _dadd(size_t n, double *a, double b) {
  size_t i;
  for (i=0; i<n; ++i) a[i] += b;
}

// a(:,*) = a(:,*) .+ b
void _dgcmv(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa, *pb;
  for (i=0,pa=a; i<n; ++i)
    for (j=0,pb=b; j<m; ++j, ++pa, ++pb)
      *pa += *pb;
}

// a(*,:) = a(*,:) .+ b
void _dgrmv(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa =a, *pb =b;
  for (i=0; i<n; ++i,++pb)
    for (j=0; j<m; ++j, ++pa)
      *pa += *pb;
}

// a = diag(b) * a
void _dgcms(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa = a, *pb;
  for (i=0; i<n; ++i)
    for (j=0, pb=b; j<m; ++j, ++pa, ++pb)
      *pa *= *pb;
}

// a = a * diag(b)
void _dgrms(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa = a, *pb = b;
  for (i=0; i<n; ++i,++pb)
    for (j=0; j<m; ++j, ++pa)
      *pa *= *pb;
}

// a = diag(1./b) * a
void _dicms(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa, *pb;
  for (j=0; j<m; ++j) assert(b[j] > 0);
  for (i=0,pa=a; i<n; ++i)
    for (j=0,pb=b; j<m; ++j, ++pa, ++pb)
      *pa /= *pb;
}

// a = a * diag(1./b)
void _dirms(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa, *pb;
  for (i=0; i<n; ++i) assert(b[i] > 0);
  for (i=0,pa=a,pb=b; i<n; ++i,++pb) {
    for (j=0; j<m; ++j, ++pa)
      *pa /= *pb;
  }
}

// b(*) = sum(a(:,*))
void _dcsum(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa, *pb;
  for (i=0,pa=a,pb=b; i<n; ++i, ++pb) {
    *pb = 0;
    for (j=0; j<m; ++j, ++pa)
      *pb += *pa;
  }
}

// b(*) += sum(a(:,*))
void _dcsum2(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa, *pb;
  for (i=0,pa=a,pb=b; i<n; ++i, ++pb) {
    for (j=0; j<m; ++j, ++pa)
      *pb += *pa;
  }
}


// b(*) = sum(a(*,:))
void _drsum(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa, *pb;
  for (j=0,pb=b; j<m; ++j, ++pb)
    *pb = 0;
  for (i=0,pa=a; i<n; ++i)
    for (j=0,pb=b; j<m; ++j, ++pa, ++pb)
      *pb += *pa;
}

// b(*) += sum(a(*,:))
void _drsum2(size_t m, size_t n, double *a, double *b) {
  size_t i,j;
  double *pa, *pb;
  for (i=0,pa=a; i<n; ++i)
    for (j=0,pb=b; j<m; ++j, ++pa, ++pb)
      *pb += *pa;
}

// normalize by column
void _dcnorm(size_t m, size_t n, double *a, double *sa) {
  size_t i, j;
  double *pa;
  bool isAllocated = true;
  if (!sa) {
    isAllocated = false;
    sa = _MALLOC_DOUBLE(n);
  }
  _dcsum(m, n, a, sa);
  for (i=0; i<n; ++i) assert(sa[i] > 0);
  for (i=0, pa=sa; i<n; ++i, ++pa) {
    for (j=0; j<m; ++j, ++a) (*a) /= *pa;
  }
  if (!isAllocated) _FREE(sa);
}

// normalize by row
void _drnorm(size_t m, size_t n, double *a, double *sa) {
  size_t i, j;
  double *pa;
  bool isAllocated = true;
  if (!sa) {
    isAllocated = false;
    sa = _MALLOC_DOUBLE(m);
  }
  _drsum(m, n, a, sa);
  for (i=0; i<m; ++i) assert(sa[i] > 0);
  for (i=0; i<n; ++i) {
    pa = sa;
    for (j=0; j<m; ++j, ++a, ++pa) (*a) /= *pa;
  }
  if (!isAllocated) _FREE(sa);
}

// center by column
void _dccenter(size_t m, size_t n, double *a, double *sa) {
  size_t i, j;
  double *pa;
  bool isAllocated = true;
  if (!sa) {
    isAllocated = false;
    sa = _MALLOC_DOUBLE(n);
  }
  _dcsum(m, n, a, sa);
  cblas_dscal(n, 1./m, sa, 1);
  for (i=0, pa=sa; i<n; ++i, ++pa) {
    for (j=0; j<m; ++j, ++a) (*a) -= *pa;
  }
  if (!isAllocated) _FREE(sa);
}

// center by row
void _drcenter(size_t m, size_t n, double *a, double *sa) {
  size_t i, j;
  double *pa;
  bool isAllocated = true;
  if (!sa) {
    isAllocated = false;
    sa = _MALLOC_DOUBLE(m);
  }
  _drsum(m, n, a, sa);
  cblas_dscal(m, 1./n, sa, 1);
  for (i=0; i<n; ++i) {
    pa = sa;
    for (j=0; j<m; ++j, ++a, ++pa) (*a) -= *pa;
  }
  if (!isAllocated) _FREE(sa);
}

// c = a.*b
void _dvmul(size_t n, double *a, double *b, double *c) {
  size_t i;
  for (i=0; i<n; ++i, ++c, ++a, ++b)
    *c = (*a) * (*b);
}


void _dpdist2(int d, size_t n, size_t m, double * A, double * B, double *C) {
  size_t i, j, ki, kj; int k;
  assert(d>0 && n>0 && m>0);

  for (i=0; i<m*n; ++i) C[i] = 0;
  for (i=0; i<m; ++i)
    for (j=0; j<n; ++j)
      for (k=0, kj=j*d, ki=i*d; k<d; ++k, ++kj, ++ki)
        C[i*n + j] += (A[kj] -  B[ki]) * (A[kj] -  B[ki]);
}

void _dpdist2_sym(int d, size_t n, size_t m, double *A, int *Bi, double *C, const double *vocab) {
  size_t i, j, ki, kj; int k;
  for (i=0; i<m*n; ++i) C[i] = 0;
  for (i=0; i<m; ++i)
    for (j=0; j<n; ++j)
      for (k=0, kj=j*d, ki=Bi[i]*d; k<d; ++k, ++kj, ++ki)
        if (Bi[i] < 0)
          C[i*n + j] += A[kj]*A[kj];
        else
          C[i*n + j] += (A[kj] - vocab[ki]) * (A[kj] - vocab[ki]);
}

void _dpdist2_submat(size_t m, int *Bi, double *C,
                     const int vocab_size, const double *dist_mat) {
  size_t i; int j;
  assert(m>0);

  for (i=0; i<m; ++i)
    for (j=0; j<vocab_size; ++j)
      C[i*vocab_size + j] = dist_mat[Bi[i]*vocab_size + j];
}

void _dpdist_symbolic(int d, size_t n, size_t m, int * A, int * B, double *C,
                      const int vocab_size, const double* dist_mat) {
  size_t i,j; int k;
  assert(d>0 && n>0 && m>0);

  for (i=0; i<m*n; ++i) C[i] = 0;
  for (i=0; i<m; ++i)
    for (j=0; j<n; ++j)
      for (k=0; k<d; ++k)
        C[i*n+j] += dist_mat[A[j*d + k]*vocab_size + B[i*d + k]];
}

// inplace a -> exp(a)
void _dexp(size_t n, double *a) {
  size_t i;
  for (i=0; i<n; ++i, ++a) *a = exp(*a);
}

