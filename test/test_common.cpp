//
//  test_common.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/21/17.
//  Copyright © 2017 cyk. All rights reserved.
//
#include <glog/logging.h>
#include "gtest/gtest.h"
#include "test_common.h"


#if defined __APPLE__
std::string root_path = "/Users/yzc147/Dropbox/GMMHMM/code/hmm_aw";
#elif defined __unix__
/*
 *std::string root_path = "/storage/home/yzc147/work/hmm_aw";
 */
std::string root_path = "/home/ustc/work/hmm_aw";
#endif

void expect_same_hmm(HmmModel& md1, HmmModel& md2, double tol){
  int dim = md1.dim, numst = md1.numst;
  int dim2 = md2.dim, numst2 = md2.numst;
  EXPECT_EQ(dim, dim2);
  EXPECT_EQ(numst, numst2);
  for (int i=0; i<numst; i++) {
    EXPECT_NEAR(md1.a00[i], md2.a00[i], tol);
  }
  // compare transition matrix:
  for (int i=0; i<numst; i++) {
    for (int j=0; j<numst; j++)
      EXPECT_NEAR(md1.a[i][j], md2.a[i][j], tol);
  }
  // compare means:
  for (int i=0; i<numst; i++) {
    for (int j=0; j<dim; j++)
      EXPECT_NEAR(md1.stpdf[i].mean[j], md2.stpdf[i].mean[j], tol);
  }
  // compare sigmas:
  for (int i=0; i<numst; i++) {
    for (int m=0; m<md1.stpdf[i].dim; m++) {
      for (int n=0; n<md1.stpdf[i].dim; n++)
        EXPECT_NEAR(md1.stpdf[i].sigma[m][n], md2.stpdf[i].sigma[m][n], tol);
    }
  }
}

