//
//  mosek_solver_test.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/8/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include "gtest/gtest.h"
#include "glog/logging.h"
#include "dist_utils.h"
#include "mosek_solver.h"
#include "utils.h"


TEST(MatchByDistmatTest, CalcDistmat){
  const char filename[] = "/Users/yzc147/Dropbox/GMMHMM/code/hmmaw/data/test/hmm2.in";
  HmmModel *hmm1;
  hmm1=(HmmModel *)calloc(1,sizeof(HmmModel));
  hmm_read(hmm1, filename);
  int n = hmm1->numst, m = hmm1->numst;
  double *C = (double *)calloc(n*m, sizeof(double));
  calc_distmat(hmm1, hmm1, C);
  double gt_dist[] = {0.0, 0.02, 0.08,\
    0.02, 0.0, 0.02,\
    0.08, 0.02, 0.0};
  double* x = (double*) calloc(n*m, sizeof(double));
  match_by_distmat(n, m, C, hmm1->a00, hmm1->a00, x, NULL);
  print_mat(x, n, m);
}
