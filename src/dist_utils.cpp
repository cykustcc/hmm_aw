//
//  dist_utils.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/3/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include "dist_utils.h"
#include "math.h"
#include <assert.h>

// Assume GaussModel A and B has diagnal sigmas.
// C[j][i] = C[j*m + i], C is an n by m matrix. Row Major.
void pdist2_hmm(size_t d, size_t n, size_t m, HmmModel& A, HmmModel& B,
                double *C){
  assert(d>0 && n>0 && m>0);
  for (size_t i=0; i<m*n; ++i) C[i] = 0;
  for (size_t i=0; i<m; ++i)
    for (size_t j=0; j<n; ++j)
      for (size_t k=0; k<d; ++k){
      C[j*m + i] += (A.stpdf[i].mean[k] - B.stpdf[j].mean[k]) * (A.stpdf[i].mean[k] - B.stpdf[j].mean[k]);
      C[j*m + i] += A.stpdf[i].sigma[k][k] + B.stpdf[j].sigma[k][k] - 2*sqrt(A.stpdf[i].sigma[k][k]*B.stpdf[j].sigma[k][k]);
      }
}


void calc_distmat(HmmModel& hmm1, HmmModel& hmm2, double *C){
  assert(hmm1.dim == hmm2.dim);
  size_t n = hmm1.numst, m = hmm2.numst, d = hmm1.dim;
  pdist2_hmm(d, n, m, hmm1, hmm2, C);
}
