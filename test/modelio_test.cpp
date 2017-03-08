//
//  modelio_test.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/7/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include "gtest/gtest.h"
#include "glog/logging.h"
#include "hmm.h"
#include <utils/blas_utils.h>

TEST(ModelIoTest, HmmRead){
  const char filename[] = "/Users/yzc147/Dropbox/GMMHMM/code/hmmaw/data/test/hmm1.in";
  HmmModel *hmm1;
  hmm1=(HmmModel *)calloc(1,sizeof(HmmModel));
  hmm_read(hmm1, filename);
  EXPECT_EQ(hmm1->dim, 2);
  EXPECT_EQ(hmm1->numst, 3);
  // Test inverse matrix function.
  for(int i=0; i<hmm1->numst; i++){
    for (int j=0; j<hmm1->dim; j++) {
      double dot_prod = cblas_ddot(hmm1->dim, hmm1->stpdf[i]->sigma[j], 1, hmm1->stpdf[i]->sigma_inv[j], 1);
      EXPECT_NEAR(dot_prod, 1.0, 0.01);
    }
  }
  print_model(hmm1, stdout);
  freehmm(&hmm1);
}

TEST(ModelIoTest, HmmWrite){
  const char filenamein[] = "/Users/yzc147/Dropbox/GMMHMM/code/hmmaw/data/test/hmm1.in";
  const char filenameout[] = "/Users/yzc147/Dropbox/GMMHMM/code/hmmaw/data/test/hmm1.out";
  HmmModel *hmm1;
  hmm1=(HmmModel *)calloc(1,sizeof(HmmModel));
  hmm_read(hmm1, filenamein);
  EXPECT_EQ(hmm1->dim, 2);
  EXPECT_EQ(hmm1->numst, 3);
  hmm_write(hmm1, filenameout);
  freehmm(&hmm1);
}

