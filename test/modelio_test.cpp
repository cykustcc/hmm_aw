//
//  modelio_test.cpp
//  hmmjia
//
//  Created by Yukun Chen on 3/7/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include "gtest/gtest.h"
#include "glog/logging.h"
#include "hmm.h"
#include <utils/blas_utils.h>

TEST(ModelIoTest, HmmRead){
  const char filename[] = "/Users/yzc147/Dropbox/GMMHMM/code/hmmjia/hmmjia/test/hmm1.in";
  HmmModel *hmm1;
  hmm1=(HmmModel *)calloc(1,sizeof(HmmModel));
  hmm_read(hmm1, filename);
  EXPECT_EQ(hmm1->dim, 2);
  EXPECT_EQ(hmm1->numst, 3);
  print_model(hmm1, stdout);
  freehmm(&hmm1);
}
