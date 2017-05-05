//
//  hmmdist_test.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/27/17.
//  Copyright © 2017 cyk. All rights reserved.
//
#include "gtest/gtest.h"
#include "glog/logging.h"
#include <hmm.h>
#include "test_common.h"

class HMMDistTest: public ::testing::Test{
protected:
  HmmModel gt_hmm;
  virtual void SetUp(){
    std::string filenamein = root_path + "/data/test/hmm_forhmmfittest.in";
    gt_hmm.read_model(filenamein);
  }
  virtual void TearDown()
  {
    
  }
};

TEST_F(HMMDistTest, GenSeq){
  int n = 800;
  std::vector<std::vector<float>> seq(1, std::vector<float>(n*gt_hmm.dim, 0.0));
  gt_hmm.gen_seq(seq[0], n, false);
//  for (int i=0; i<n*gt_hmm.dim; i++) {
//    std::cout<<seq[0][i]<<" ";
//    if (i%gt_hmm.dim == gt_hmm.dim-1) {
//      std::cout<<std::endl;
//    }
//  }
  int dim = gt_hmm.dim, numst = gt_hmm.numst;
  std::vector<int> stcls;
  HmmModel est_hmm;
  est_hmm.resize(dim, numst, numst, stcls);
  std::vector<int> len = {n};
  std::vector<double> loglikehd = {0.0};
  double lhsumpt;
  std::vector<double> tmpwt;
  est_hmm.hmmfit(seq, 1, len, loglikehd, lhsumpt, EPSILON, tmpwt, false);
  est_hmm.print_model("");
  
  expect_same_hmm(gt_hmm, est_hmm, 0.4);
}

TEST_F(HMMDistTest, DistKL){
  double dist = gt_hmm.dist_KL(gt_hmm, 200);
  EXPECT_NEAR(dist, 0.0, 0.01);
}

TEST_F(HMMDistTest, DistMAW){
  double dist = gt_hmm.dist_MAW(gt_hmm, 0.2);
  EXPECT_NEAR(dist, 0.0, 0.01);
}

TEST_F(HMMDistTest, DistIAW){
  double dist = gt_hmm.dist_IAW(gt_hmm, 0.2, 100, true);
  EXPECT_NEAR(dist, 0.0, 0.01);
}

