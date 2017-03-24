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

void same_hmm(HmmModel& md1, HmmModel& md2, double tol){
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

class HMMFitTest: public ::testing::Test{
protected:
  HmmModel est_hmm_diag, est_hmm_notdiag;
  HmmModel gt_hmm;
  virtual void SetUp(){
    std::string datafile = root_path + "data/test/seqhmmdim2st2.txt";
    int len = 200, dim = 2, numst = 2;
    std::vector<double> loglikehd(1, 0.0);
    
    // TODO: add a function to read the sequence data!
    std::vector<std::vector<float>> u(1, std::vector<float>(dim*len, 0.0)); //since nseq = 1
    FILE *fp = fopen(datafile.c_str(), "r");
    for (int i=0; i<len; i++) {
      fscanf(fp, "%f, %f", &u[0][2*i], &u[0][2*i+1]);
    }
    double lhsum;
    std::vector<double> tmpwt;
    std::vector<int> stcls;
    int numcls = stcls.size() > 0 ? *std::max_element(stcls.begin(), stcls.end()) : numst;
    std::vector<int> lens; lens.push_back(len);
    est_hmm_notdiag.resize(dim, numst, numcls, stcls);
    est_hmm_diag.resize(dim, numst, numcls, stcls);
    hmmfit(est_hmm_notdiag, u, 1, lens, loglikehd, lhsum, EPSILON, tmpwt, false);
    hmmfit(est_hmm_diag, u, 1, lens, loglikehd, lhsum, EPSILON, tmpwt, true);
    
    std::string filenamein = root_path + "/data/test/hmm_forhmmfittest.in";
    gt_hmm.read_model(filenamein);
  }
  virtual void TearDown()
  {
    
  }
};

TEST_F(HMMFitTest, NotForceSigmaBeDiagnal){
  same_hmm(gt_hmm, est_hmm_notdiag, 0.4);
  est_hmm_notdiag.print_model("");
}

TEST_F(HMMFitTest, ForceSigmaBeDiagnal){
  same_hmm(gt_hmm, est_hmm_diag, 0.4);
  est_hmm_diag.print_model("");
}
