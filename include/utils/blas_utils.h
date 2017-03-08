//
//  blas_utils.h
//  hmmjia
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

#endif /* blas_utils_h */
