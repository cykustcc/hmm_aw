//
//  estimate_test.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/11/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include "gtest/gtest.h"
#include "glog/logging.h"
#include "test_common.h"
#include "dist_utils.h"
#include "mosek_solver.h"
#include "hmm.h"
#include <cstring>

void same_hmm(HmmModel * md1, HmmModel * md2, double tol){
  int dim = md1->dim, numst = md1->numst;
  int dim2 = md2->dim, numst2 = md2->numst;
  EXPECT_EQ(dim, dim2);
  EXPECT_EQ(numst, numst2);
  for (int i=0; i<numst; i++) {
    EXPECT_NEAR(md1->a00[i], md2->a00[i], tol);
  }
  // compare transition matrix:
  for (int i=0; i<numst; i++) {
    for (int j=0; j<numst; j++)
      EXPECT_NEAR(md1->a[i][j], md2->a[i][j], tol);
  }
  // compare means:
  for (int i=0; i<numst; i++) {
    for (int j=0; j<dim; j++)
      EXPECT_NEAR(md1->stpdf[i]->mean[j], md2->stpdf[i]->mean[j], tol);
  }
  // compare sigmas:
  for (int i=0; i<numst; i++) {
    for (int m=0; m<md1->stpdf[i]->dim; m++) {
      for (int n=0; n<md1->stpdf[i]->dim; n++)
        EXPECT_NEAR(md1->stpdf[i]->sigma[m][n], md2->stpdf[i]->sigma[m][n], tol);
    }
  }
}

class HMMFitTest: public ::testing::Test{
protected:
  HmmModel *est_hmm_diag, *est_hmm_notdiag;
  HmmModel *gt_hmm;
  virtual void SetUp(){
    char datafile[] = "data/test/seqhmmdim2st2.txt";
    char filename[40];
    std::strcpy(filename, root_path);
    std::strcat(filename, datafile);
    int len = 200, dim = 2, numst = 2;
    double * loglikehd=(double *)calloc(1,sizeof(double));
    
    FILE *fp = fopen(filename, "r");
    float *data = (float *)calloc(len*dim, sizeof(float));
    for (int i=0; i<len; i++) {
      fscanf(fp, "%f, %f", &data[2*i], &data[2*i+1]);
    }
    float **u= &data;
    double lhsum;
    est_hmm_notdiag=(HmmModel *)calloc(1,sizeof(HmmModel));
    hmmfit(u, 1, &len, dim, est_hmm_notdiag, numst, NULL, loglikehd, &lhsum, EPSILON, NULL, false);
    est_hmm_diag=(HmmModel *)calloc(1,sizeof(HmmModel));
    hmmfit(u, 1, &len, dim, est_hmm_diag, numst, NULL, loglikehd, &lhsum, EPSILON, NULL, true);
    
    const char filenamein[] = "/Users/yzc147/Dropbox/GMMHMM/code/hmmaw/data/test/hmm_forhmmfittest.in";
    gt_hmm=(HmmModel *)calloc(1,sizeof(HmmModel));
    hmm_read(gt_hmm, filenamein);
  }
  virtual void TearDown()
  {
    
  }
};

TEST_F(HMMFitTest, NotForceSigmaBeDiagnal){
  same_hmm(gt_hmm, est_hmm_notdiag, 0.4);
  print_model(est_hmm_notdiag, stdout);
}

TEST_F(HMMFitTest, ForceSigmaBeDiagnal){
  same_hmm(gt_hmm, est_hmm_diag, 0.4);
  print_model(est_hmm_diag, stdout);
}
