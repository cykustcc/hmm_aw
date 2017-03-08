//
//  dist_utils.cpp
//  hmmjia
//
//  Created by Yukun Chen on 3/3/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include "dist_utils.h"
#include "math.h"
#include <assert.h>

// Assume GaussModel A and B has diagnal sigmas.
void pdist2_hmm(size_t d, size_t n, size_t m, HmmModel* A, HmmModel* B, double *C) {
  size_t i, j; int k;
  assert(d>0 && n>0 && m>0);
  for (i=0; i<m*n; ++i) C[i] = 0;
  for (i=0; i<m; ++i)
  for (j=0; j<n; ++j)
    for (k=0; k<d; ++k){
      C[i*n + j] += (A->stpdf[i]->mean[k] - B->stpdf[j]->mean[k]) * (A->stpdf[i]->mean[k] - B->stpdf[j]->mean[k]);
      C[i*n + j] += A->stpdf[i]->sigma[k][k] + B->stpdf[j]->sigma[k][k] - 2*sqrt(A->stpdf[i]->sigma[k][k]*B->stpdf[j]->sigma[k][k]);
    }
}


void calc_distmat(HmmModel* hmm1, HmmModel* hmm2, double *C){
  assert(hmm1->dim == hmm2->dim);
  size_t n = hmm1->numst, m = hmm2->numst, d = hmm1->dim;
  pdist2_hmm(d, n, m, hmm1, hmm2, C);
}
