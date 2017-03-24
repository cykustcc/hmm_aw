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
  std::string filename("/Users/yzc147/Dropbox/GMMHMM/code/hmmaw/data/test/hmm1.in");
  HmmModel hmm1;
  hmm1.read_model(filename);
  EXPECT_EQ(hmm1.dim, 2);
  EXPECT_EQ(hmm1.numst, 3);
  // Test inverse matrix function.
  for(int i=0; i<hmm1.numst; i++){
    for (int j=0; j<hmm1.dim; j++) {
      double dot_prod = cblas_ddot(hmm1.dim,
                                   (double *) &hmm1.stpdf[i].sigma[j],
                                   1,
                                   (double *) &hmm1.stpdf[i].sigma_inv[j],
                                   1);
      EXPECT_NEAR(dot_prod, 1.0, 0.01);
    }
  }
  hmm1.print_model("");
}

TEST(ModelIoTest, HmmWrite){
  std::string filenamein("/Users/yzc147/Dropbox/GMMHMM/code/hmmaw/data/test/hmm1.in");
  std::string filenameout("/Users/yzc147/Dropbox/GMMHMM/code/hmmaw/data/test/hmm1.out");
  HmmModel hmm1;
  hmm1.read_model(filenamein);
  EXPECT_EQ(hmm1.dim, 2);
  EXPECT_EQ(hmm1.numst, 3);
  hmm1.write_model(filenameout);
}

