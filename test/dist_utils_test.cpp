//
//  dist_utils_test.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/7/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include "gtest/gtest.h"
#include "glog/logging.h"
#include "dist_utils.h"
#include "test_common.h"


TEST(DistUtilTest, CalcDistmat){
  std::string filename(root_path + "/data/test/hmm2.in");
  HmmModel hmm1;
  hmm1.read_model(filename);
  double *C = (double *)calloc(hmm1.numst*hmm1.numst, sizeof(double));
  calc_distmat(hmm1, hmm1, C);
  double gt_dist[] = {0.0, 0.02, 0.08,\
                     0.02, 0.0, 0.02,\
                     0.08, 0.02, 0.0};
  for (int i=0; i<hmm1.numst*hmm1.numst; i++) {
    EXPECT_NEAR(C[i], gt_dist[i], 0.0001);
//    LOG(INFO)<<C[i];
  }
}
