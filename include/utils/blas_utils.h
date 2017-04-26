//
//  blas_utils.h
//  hmm_aw
//
//  modified from https://github.com/bobye/d2_kmeans Copyright @ Jianbo Ye
//

#ifndef blas_utils_h
#define blas_utils_h

#ifdef __BLAS_LEGACY__
#include <math.h>
#include "utils/cblas.h"
#define _MALLOC_DOUBLE(x)       (double *) malloc( (x) *sizeof(double))
#define _MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _MALLOC_SIZE_T(x)       (size_t *) malloc( (x) *sizeof(size_t))
#define _CALLOC_DOUBLE(x)       (double *) calloc( (x) , sizeof(double))
#define _CALLOC_INT(x)       (int *) calloc( (x) , sizeof(int))
#define _CALLOC_SIZE_T(x)       (size_t *) calloc( (x) , sizeof(size_t))
#define _FREE(x)         free(x)

#elif defined __APPLE__
#include <Accelerate/Accelerate.h>
#define _MALLOC_DOUBLE(x)       (double *) malloc( (x) *sizeof(double))
#define _MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _MALLOC_SIZE_T(x)       (size_t *) malloc( (x) *sizeof(size_t))
#define _CALLOC_DOUBLE(x)       (double *) calloc( (x) , sizeof(double))
#define _CALLOC_INT(x)       (int *) calloc( (x) , sizeof(int))
#define _CALLOC_SIZE_T(x)       (size_t *) calloc( (x) , sizeof(size_t))
#define _FREE(x)         free(x)

#elif defined __USE_MKL__
#include <mkl.h>
#define _MALLOC_DOUBLE(x)       (double *) mkl_malloc( (x) *sizeof(double), 16)
#define _MALLOC_INT(x)       (int *) mkl_malloc( (x) *sizeof(int), 16)
#define _MALLOC_SIZE_T(x)       (size_t *) mkl_malloc( (x) *sizeof(size_t), 16)
#define _CALLOC_DOUBLE(x)       (double *) mkl_calloc( (x) , sizeof(double), 16)
#define _CALLOC_INT(x)       (int *) mkl_calloc( (x) , sizeof(int), 16)
#define _CALLOC_SIZE_T(x)       (size_t *) mkl_calloc( (x) , sizeof(size_t), 16)
#define _FREE(x)         mkl_free(x)
#endif


template <typename DType>
void Gemm(const CBLAS_TRANSPOSE TransA,
          const CBLAS_TRANSPOSE TransB,
          const int M, /* # of rows of A and C */
          const int N, /* # of cols of A and rows of B */
          const int K, /* # of cols of B and C */
          const DType alpha,
          const DType* A,
          const DType* B,
          const DType beta,
          DType* C,
          bool isdouble = false /* float or double, default float*/) {
  int lda = (TransA == CblasNoTrans) ? K : M;
  int ldb = (TransB == CblasNoTrans) ? N : K;
  if(!isdouble)
    cblas_sgemm(CblasRowMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb,
              beta, C, N);
  else
    cblas_dgemm(CblasRowMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb,
                beta, C, N);
}

#endif /* blas_utils_h */
