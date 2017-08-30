//
//  optimaltransport.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 4/28/17.
//  Copyright © 2017 cyk. All rights reserved.
//
#include "stdlib.h"
#include <cassert>

#include <optimaltransport.h>
#include "utils/blas_utils.h"
#include "utils/blas_like.h"


#define ROUNDOFF (1E-9)
#define BADMM_TOL (1E-4)

double match_by_distmat_BADMM(int n,
                              int m,
                              double *C,
                              double* wX,
                              double* wY,
                              /** OUT **/ double *x,
                              /** OUT **/ double *lambda,
                              int rho,
                              int max_iter){
  assert(rho > 0);
  //  Pi2=W1*W2';
  //  Lambda=zeros(size(Pi2));
  double *Pi1 = x;
  double *Pi2 = _MALLOC_DOUBLE(n*m);
  double *Lambda = _MALLOC_DOUBLE(n*m);
  //  C=C/rho;
  //  xi=exp(-C);
  double *C_div_rho = _MALLOC_DOUBLE(n*m);
  double *xi = _MALLOC_DOUBLE(n*m);
  cblas_dcopy(n*m, C, 1, C_div_rho, 1);
  cblas_dscal(n*m, 1.0/rho, C_div_rho, 1);
  for (int i=0; i<n*m; ++i, ++xi) *xi = exp(-*C_div_rho);
  double *tmp = _MALLOC_DOUBLE(n*m);
  for (int iter=0; iter<max_iter; iter++) {
    // tmp = exp(Lambda);
    for (int j=0; j<n*m; ++j, ++xi) *tmp = exp(*Lambda);
    // Pi1=Pi2 .* xi ./tmp  + eps;
    for (int j=0; j<n*m; j++) Pi1[j] = Pi2[j] * xi[j] / tmp[j] + ROUNDOFF;
    // Pi1=bsxfun(@times, Pi1, W1 ./sum(Pi1, 2));
    _dcnorm(n, m, Pi1, NULL);
    _dgcmv(n, m, Pi1, wX);
    // Pi2=Pi1 .* tmp + eps;
    for (int j=0; j<n*m; j++) Pi2[j] = Pi1[j] * tmp[j] + ROUNDOFF;
    // Pi2=bsxfun(@times, Pi2', W2 ./sum(Pi2, 1)')'; % memory overheads in matrix transpose
    _drnorm(n, m, Pi2, NULL);
    _dgrmv(n, m, Pi1, wY);
    // Lambda=Lambda + Pi1 - Pi2;
    for (int j=0; j<n*m; ++j, ++xi) Lambda[j] += Pi1[j] - Pi2[j];
    // if mod(i,100)==0 && max(abs(Pi1(:)-Pi2(:)))<eps
    //   break
    // end
    if ((iter%100==99 ))  {
      double maxdiff = 0.0;
      for (int j=0; j<n*m; ++j, ++xi)
        maxdiff =std::max(maxdiff, std::abs(Pi1[j] - Pi2[j]));
      if (maxdiff < BADMM_TOL) {
        break;
      }
    }
  }
  
  _FREE(Pi1);
  _FREE(Pi2);
  _FREE(Lambda);
  _FREE(C_div_rho);
  _FREE(xi);
  return 0.0;
}
