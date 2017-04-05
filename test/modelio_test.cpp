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
#include "test_common.h"
#include <utils/blas_utils.h>

TEST(ModelIoTest, HmmRead){
  std::string filename(root_path + "/data/test/hmm1.in");
  HmmModel hmm1;
  hmm1.read_model(filename);
  EXPECT_EQ(hmm1.dim, 2);
  EXPECT_EQ(hmm1.numst, 3);
  // Test inverse matrix function.
  for(int i=0; i<hmm1.numst; i++){
    for (int j=0; j<hmm1.dim; j++) {
      double dot_prod = 0;
      for (int k=0; k<hmm1.dim; k++) {
        dot_prod += hmm1.stpdf[i].sigma[j][k] * hmm1.stpdf[i].sigma_inv[j][k];
      }
      EXPECT_NEAR(dot_prod, 1.0, 0.01);
    }
  }
  hmm1.print_model("");
}

TEST(ModelIoTest, HmmWrite){
  std::string filenamein(root_path + "/data/test/hmm1.in");
  std::string filenameout(root_path + "/data/test/hmm1.out");
  HmmModel hmm1;
  hmm1.read_model(filenamein);
  EXPECT_EQ(hmm1.dim, 2);
  EXPECT_EQ(hmm1.numst, 3);
  hmm1.write_model(filenameout);
}

