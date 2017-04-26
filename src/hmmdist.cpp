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
#include "matrix.h"
#include "dist_utils.h"
#include "mosek_solver.h"


double HmmModel::dist_KL(HmmModel &hmm2, int sample_size, bool diag){
  double res = 0.0;
  std::vector<std::vector<float>> seq(1, std::vector<float>(dim*sample_size, 0.0));
  gen_seq(seq[0], sample_size, diag);
  std::vector<int> len = {sample_size};
  std::vector<double> wt;
  res = comploglike(seq, 1, len, wt);
  res -= hmm2.comploglike(seq, 1, len, wt);
  res /= sample_size;
  return res;
}


void normalize_row(double* x, double* x_r, int m, int n){
  double *row_sums = (double *)malloc(m*sizeof(double));
  for (int i=0; i<m; i++) {
    row_sums[i] = 0;
    for (int j=0; j<n; j++) {
      row_sums[i] += x[i*n+j];
    }
    if (fabs(row_sums[i] - 0.0) < 1e-5) {
      row_sums[i] = 1.0;
    }
  }
  
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      x_r[i*n+j] = x[i*n+j] / row_sums[i];
    }
  }
}

void normalize_col(double* x, double* x_c, int m, int n){
  double *col_sums = (double *)malloc(m*sizeof(double));
  for (int i=0; i<n; i++) {
    col_sums[i] = 0;
    for (int j=0; j<m; j++) {
      col_sums[i] += x[i+j*m];
    }
    if (fabs(col_sums[i] - 0.0) < 1e-5) {
      col_sums[i] = 1.0;
    }
  }
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      x_c[i+j*m] = x[i*n+j] / col_sums[i];
    }
  }
}

double HmmModel::dist_transmat_MAW(HmmModel &hmm2, double* C, double* x){
  double *x_w = (double *)malloc(numst*hmm2.numst*sizeof(double));
  double *x_r = (double *)malloc(numst*hmm2.numst*sizeof(double));
  normalize_row(x, x_r, numst, hmm2.numst);
  normalize_col(x, x_r, numst, hmm2.numst);
  free(x_w);
  free(x_r);
  return 0.0;
}


double HmmModel::dist_MAW(HmmModel &hmm2){
  double *C = (double *)malloc(numst*hmm2.numst*sizeof(double));
  double *x = (double *)malloc(numst*hmm2.numst*sizeof(double));
  calc_distmat(*this, hmm2, C);
  double lambda = 1;
  
  solver_setup();
  double val = match_by_distmat(numst,
                                hmm2.numst,
                                C,
                                a00,
                                hmm2.a00,
                                x,
                                &lambda);
  solver_release();
  
  free(C);
  free(x);
  return val;
}

double HmmModel::dist_IAW(HmmModel &hmm2){
  return 0.0;
}

// based on current state('s gaussian), generate the idx th point of the sequence.
// i.e. seq[idx*dim:idx*(dim+1)-1] are filled in.
void HmmModel::gauss_sample(std::vector<float> &seq, /*output*/
                            int idx,
                            int state,
                            std::vector<std::vector<double>> &mt_S,
                            bool diag){
  std::vector<double> res(dim, 0.0);
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      res[i] += mt_S[i][j]*seq[idx*dim+j];
    }
    res[i] += stpdf[state].mean[i];
    seq[idx*dim+i] = res[i];
  }
}

void HmmModel::gen_seq(std::vector<float> &seq, int n, bool diag){
  assert(seq.size() == n*dim);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::map<int, int> m;
  std::vector<std::discrete_distribution<>> state_dd;
  std::discrete_distribution<> init_state_dd(this->a00.begin(), this->a00.end());
  std::normal_distribution<> gauss_std(0,1);
//  std::map<int, int> hist;// to be deleted
  for (int i=0; i<numst; i++) {
    std::discrete_distribution<> tmp(this->a[i].begin(), this->a[i].end());
    state_dd.push_back(tmp);
  }
  int cur_state = init_state_dd(gen);
  std::vector<std::vector<double>> mt_S(dim, std::vector<double>(dim, 0.0));
  std::vector<std::vector<std::vector<double>>> mt_Ss;
  for (int i=0; i<numst; i++) {
    mt_Ss.push_back(mt_S);
    cholesky_decomp(stpdf[i].sigma,
                    mt_Ss[i],
                    diag);
  }
  for(int i = 0; i < n; i++){
    double cumulate = 0.0;
    cur_state = state_dd[cur_state](gen);
    for (int j = 0; j < dim; j++) {
      seq[i*dim+j] = gauss_std(gen);
//      ++hist[std::round(seq[i*dim+j])];
    }
    gauss_sample(seq, i, cur_state, mt_Ss[cur_state], diag);
    m[cur_state]++;
  }
  // to be deleted
//  for(auto p : hist) {
//    std::cout << std::fixed
//    << p.first << ' ' << std::string(p.second/10, '*') << '\n';
//  }
  for(auto p : m) {
    std::cout << p.first << " generated " << p.second << " times\n";
  }
};
