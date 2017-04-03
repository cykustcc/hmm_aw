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
  int n = 200;
  std::vector<float> seq(n, 0.0);
  gt_hmm.gen_seq(seq, n);
}

