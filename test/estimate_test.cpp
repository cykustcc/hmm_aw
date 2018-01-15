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
#include "hmm.h"
#include <cstring>
#include "matrix.h"


class HMMFitTest: public ::testing::Test{
protected:
  HmmModel est_hmm_diag, est_hmm_notdiag;
  HmmModel gt_hmm;
  virtual void SetUp(){
    std::string datafile = root_path + "/data/test/seqhmmdim2st2.txt";
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
    est_hmm_notdiag.hmmfit(u, 1, lens, loglikehd, lhsum, EPSILON, tmpwt, false);
    est_hmm_diag.hmmfit(u, 1, lens, loglikehd, lhsum, EPSILON, tmpwt, true);

    std::string filenamein = root_path + "/data/test/hmm_forhmmfittest.in";
    gt_hmm.read_model(filenamein);
  }
  virtual void TearDown()
  {

  }
};

TEST_F(HMMFitTest, NotForceSigmaBeDiagnal){
  expect_same_hmm(gt_hmm, est_hmm_notdiag, 0.4);
  est_hmm_notdiag.print_model("");
}

TEST_F(HMMFitTest, ForceSigmaBeDiagnal){
  expect_same_hmm(gt_hmm, est_hmm_diag, 0.4);
  est_hmm_diag.print_model("");
}

