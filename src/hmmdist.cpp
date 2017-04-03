//
//  hmmdist.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/27/17.
//  Copyright © 2017 cyk. All rights reserved.
//
#include "glog/logging.h"
#include "hmm.h"
#include <random>
#include <cassert>
#include <map>
#include "utils/lapacke.h"

double HmmModel::dist_KL(HmmModel &hmm2){
  return 0.0;
}

double HmmModel::dist_MAW(HmmModel &hmm2){
  return 0.0;
}

double HmmModel::dist_IAW(HmmModel &hmm2){
  return 0.0;
}

void HmmModel::gen_seq(std::vector<float> &seq, int n){
  assert(seq.size() == n);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::map<int, int> m;
  std::vector<std::discrete_distribution<>> state_dd;
  std::discrete_distribution<> init_state_dd(this->a00.begin(), this->a00.end());
  for (int i=0; i<numst; i++) {
    std::discrete_distribution<> tmp(this->a[i].begin(), this->a[i].end());
    state_dd.push_back(tmp);
  }
  int cur_state = init_state_dd(gen);
  for(int i = 0; i < n; i++){
    double cumulate = 0.0;
    cur_state = state_dd[cur_state](gen);
    
    m[cur_state]++;
  }
  for(auto p : m) {
    std::cout << p.first << " generated " << p.second << " times\n";
  }
};
