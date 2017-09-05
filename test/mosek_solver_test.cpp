//
//  mosek_solver_test.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/8/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include "gtest/gtest.h"
#include "glog/logging.h"

#include "test_common.h"
#include "dist_utils.h"
#include "mosek_solver.h"


TEST(MatchByDistmatTest, MatchByDistmat){
  std::string filename = root_path + "/data/test/hmm2.in";
  HmmModel hmm1;
  hmm1.read_model(filename);
  int n = hmm1.numst, m = hmm1.numst;
  double *C = (double *)calloc(n*m, sizeof(double));
  calc_distmat(hmm1, hmm1, C);
  double gt_dist[] = {0.0, 0.02, 0.08,\
    0.02, 0.0, 0.02,\
    0.08, 0.02, 0.0};
  for (int i=0; i<n*m; i++) {
    EXPECT_NEAR(C[i], gt_dist[i], 0.0001);
  }
  double* match = (double*) calloc(n * m, sizeof(double));
  double* a00 = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    a00[i] = hmm1.a00[i];
  }
  solver_setup();
  double d = match_by_distmat(n, m, C, a00, a00, match, NULL);
  solver_release();
  EXPECT_NEAR(d, 0, 0.0001);
  double gt_match[3][3] = {{0.5,0.0,0.0},
                        {0.0,0.2,0.0},
                        {0.0,0.0,0.3}};
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      EXPECT_NEAR(match[i*m+j], gt_match[i][j], 0.0001);
    }
  }
  print_matrix(match, n, m);
}
