//
//  hmmdist.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/27/17.
//  Copyright © 2017 cyk. All rights reserved.
//
#include "glog/logging.h"
#include "hmm.h"
#include <cassert>
#include <map>
#include <random>
#include <vector>
//#include "utils/lapacke.h"
#include "matrix.h"
#include "dist_utils.h"
#include "mosek_solver.h"
#include "utils/blas_utils.h"
#include "utils/blas_like.h"
#include "utils/common.h"


double HmmModel::dist_KL(HmmModel &hmm2, int sample_size, bool diag) {
  double res = 0.0;
  std::vector<std::vector<float>> seq(
      1, std::vector<float>(dim * sample_size, 0.0));
  gen_seq(seq[0], sample_size, diag);
  std::vector<int> len = {sample_size};
  std::vector<double> wt;
  res = comploglike(seq, 1, len, wt);
  res -= hmm2.comploglike(seq, 1, len, wt);
  res /= sample_size;
  return res;
}


void normalize_row(double* x, double* x_r, int m, int n) {
  double *row_sums = (double *)malloc(m*sizeof(double));
  for (int i = 0; i < m; i++) {
    row_sums[i] = 0;
    for (int j = 0; j < n; j++) {
      row_sums[i] += x[i*n+j];
    }
    if (fabs(row_sums[i] - 0.0) < 1e-5) {
      row_sums[i] = 1.0;
    }
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      x_r[i*n+j] = x[i*n+j] / row_sums[i];
    }
  }
}

void normalize_col(double* x, double* x_c, int m, int n) {
  double *col_sums = (double *)malloc(m*sizeof(double));
  for (int i = 0; i < n; i++) {
    col_sums[i] = 0;
    for (int j = 0; j < m; j++) {
      col_sums[i] += x[i+j*m];
    }
    if (fabs(col_sums[i] - 0.0) < 1e-9) {
      col_sums[i] = 1.0;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      x_c[i+j*m] = x[i*n+j] / col_sums[i];
    }
  }
}

double HmmModel::dist_transmat_MAW(HmmModel &hmm2, double* C, double* x) {
  double res = 0.0;
  int numst_1 = numst, numst_2 = hmm2.numst;
  double *x_c = (double *) malloc(numst_1 * numst_2 * sizeof(double));
  double *x_r = (double *) malloc(numst_1 * numst_2 * sizeof(double));
  // make transmat in format usable to cblas.
  double *transmat = (double *)malloc(numst_1 * numst_1 * sizeof(double));
  double *transmat2 = (double *)malloc(numst_2 * numst_2 * sizeof(double));
  get_transmat(transmat);
  hmm2.get_transmat(transmat2);

  normalize_row(x, x_r, numst_1, numst_2);
  normalize_col(x, x_c, numst_1, numst_2);

  double *transmat1_to_2 =
      (double *) malloc(numst_2 * numst_2 * sizeof(double));
  double *transmat1_to_2_tmp =
      (double *) malloc(numst_2 * numst_1 * sizeof(double));
  //transmat1_to_2_tmp = 1.0*x_c'*transmat + 0.0*transmat1_to_2_tmp
  //  (n_2 * n_1)    (n_1 * n_2)' * (n_1 * n_1)
  Gemm(CblasTrans,
       CblasNoTrans,
       numst_2, /* # of rows of x_c' and transmat1_to_2_tmp (A and C) */
       numst_1, /* # of cols of x_c' (A) and rows of transmat (B) */
       numst_1, /* # of cols of transmat and transmat1_to_2_tmp (B and C) */
       1.0,
       x_c,
       transmat,
       0.0,
       transmat1_to_2_tmp,
       true);
  //transmat1_to_2 = 1.0*transmat1_to_2_tmp*x_r + 0.0*transmat1_to_2
  //  (n_2 * n_2)    (n_2 * n_1) * (n_1 * n_2)
  Gemm(CblasNoTrans,
       CblasNoTrans,
       numst_2, /* # of rows of transmat1_to_2_tmp and x_r (A and C) */
       numst_1, /* # of cols of transmat1_to_2_tmp (A) and rows of x_r' (B) */
       numst_2, /* # of cols of x_r and transmat1_to_2 (B and C) */
       1.0,
       transmat1_to_2_tmp,
       x_r,
       0.0,
       transmat1_to_2,
       true);

  double *transmat2_to_1 =
      (double *) malloc(numst_1 * numst_1 * sizeof(double));
  double *transmat2_to_1_tmp =
      (double *) malloc(numst_1 * numst_2 * sizeof(double));
  //transmat2_to_1_tmp = 1.0*x_r*transmat2 + 0.0*transmat2_to_1_tmp
  //  (n_1 * n_2)    (n_1 * n_2) * (n_2 * n_2)
  Gemm(CblasNoTrans,
       CblasNoTrans,
       numst_1, /* # of rows of x_r and transmat1_to_2_tmp (A and C) */
       numst_2, /* # of cols of x_r (A) and rows of transmat (B) */
       numst_2, /* # of cols of transmat2 and transmat1_to_2_tmp (B and C) */
       1.0,
       x_c,
       transmat2,
       0.0,
       transmat1_to_2_tmp,
       true);
  //transmat2_to_1 = 1.0*transmat2_to_1_tmp*x_c' + 0.0*transmat2_to_1
  //  (n_1 * n_1)    (n_1 * n_2) * (n_1 * n_2)'
  Gemm(CblasNoTrans,
       CblasTrans,
       numst_1, /* # of rows of transmat2_to_1_tmp and x_c (A and C) */
       numst_2, /* # of cols of transmat2_to_1_tmp (A) and rows of x_c (B) */
       numst_1, /* # of cols of x_c' and transmat2_to_1_tmp (B and C) */
       1.0,
       transmat1_to_2_tmp,
       x_r,
       0.0,
       transmat1_to_2,
       true);

  for (int i = 0; i < numst_2; i++) {
    res += match_by_distmat(numst_2,
                            numst_2,
                            C,
                            transmat1_to_2+i*numst_2,
                            transmat2+i*numst_2,
                            x,
                            NULL);
  }
  for (int i = 0; i < numst_1; i++) {
    res += match_by_distmat(numst_1,
                            numst_1,
                            C,
                            transmat2_to_1+i*numst_1,
                            transmat2+i*numst_1,
                            x,
                            NULL);
  }
  free(transmat1_to_2_tmp);
  free(transmat1_to_2);
  free(transmat2_to_1_tmp);
  free(transmat2_to_1);
  free(x_c);
  free(x_r);
  return res;
}


double HmmModel::dist_MAW(HmmModel &hmm2, double alpha) {
  int numst_1 = numst, numst_2 = hmm2.numst;
  double *C = (double *) malloc(numst_1 * numst_2 * sizeof(double));
  double *match = (double *)malloc(numst_1 * numst_2 * sizeof(double));
  calc_distmat(*this, hmm2, C);
  double* hmm1_a00 = (double*) malloc(numst_1 * sizeof(double));
  double* hmm2_a00 = (double*) malloc(numst_2 * sizeof(double));
  for (int i = 0; i < numst_1; i++)
    hmm1_a00[i] = a00[i];
  for (int i = 0; i < numst_2; i++)
    hmm2_a00[i] = a00[i];
  solver_setup();
  double dist_gmm = match_by_distmat(numst_1,
                                     numst_2,
                                     C,
                                     hmm1_a00,
                                     hmm2_a00,
                                     match,
                                     NULL);
  double dist_trans = dist_transmat_MAW(hmm2, C, match);
  solver_release();

  double res = (1-alpha) * dist_gmm + alpha * dist_trans;

  free(C);
  free(match);
  return res;
}

double HmmModel::dist_IAW(HmmModel &hmm2, double alpha, int sample_size,
                          bool diag) {
  std::vector<std::vector<float>> seq(2,
                                      std::vector<float>(dim*sample_size, 0.0));
  gen_gmm(seq[0], sample_size, diag);
  hmm2.gen_gmm(seq[1], sample_size, diag);

//  print_matrix(seq);
//  match_by_distmat_BADMM();
  return 0.0;
}

// based on current state('s gaussian), generate the [idx, idx + N -1] th points of the sequence.
// i.e. seq[idx*dim:idx*(dim+1)-1] are filled in.
void HmmModel::gauss_sample(std::vector<float> &seq, /*output*/
                            int idx,
                            int N,
                            int state,
                            std::vector<std::vector<double>> &mt_S,
                            bool diag) {
  std::vector<double> res(dim, 0.0);
  for (int k = 0; k < N; k++) {
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        res[i] += mt_S[i][j] * seq[(idx + k)*dim+j];
      }
      res[i] += stpdf[state].mean[i];
      seq[(idx + k) * dim+i] = res[i];
    }
    std::fill(res.begin(), res.end(), 0.0);
  }
}

void HmmModel::gen_seq(std::vector<float> &seq, int n, bool diag) {
  assert(seq.size() == n*dim);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::map<int, int> m;
  std::vector<std::discrete_distribution<>> state_dd;
  std::discrete_distribution<> init_state_dd(this->a00.begin(),
                                             this->a00.end());
  std::normal_distribution<> gauss_std(0, 1);
//  std::map<int, int> hist;// to be deleted
  for (int i = 0; i < numst; i++) {
    std::discrete_distribution<> tmp(this->a[i].begin(), this->a[i].end());
    state_dd.push_back(tmp);
  }
  int cur_state = init_state_dd(gen);
  
  //mt_Ss stores the cholesky_decomp of sigmas.
  std::vector<std::vector<double>> mt_S(dim, std::vector<double>(dim, 0.0));
  std::vector<std::vector<std::vector<double>>> mt_Ss;
  for (int i = 0; i < numst; i++) {
    mt_Ss.push_back(mt_S);
    cholesky_decomp(stpdf[i].sigma,
                    mt_Ss[i],
                    diag);
  }

  for (int i = 0; i < n; i++) {
    cur_state = state_dd[cur_state](gen);
    for (int j = 0; j < dim; j++) {
      seq[i*dim+j] = gauss_std(gen);
//      ++hist[std::round(seq[i*dim+j])];
    }
    gauss_sample(seq, i, 1, cur_state, mt_Ss[cur_state], diag);
    m[cur_state]++;
  }
  // to be deleted
//  for(auto p : hist) {
//    std::cout << std::fixed
//    << p.first << ' ' << std::string(p.second/10, '*') << '\n';
//  }
  for (auto p : m) {
    std::cout << p.first << " generated " << p.second << " times\n";
  }
}

void HmmModel::gen_gmm(std::vector<float> &seq, int n, bool diag) {
  assert(seq.size() == n*dim);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::map<int, int> m;
  std::vector<std::discrete_distribution<>> state_dd;
  std::normal_distribution<> gauss_std(0, 1);

  std::vector<std::vector<double>> mt_S(dim, std::vector<double>(dim, 0.0));
  std::vector<std::vector<std::vector<double>>> mt_Ss;
  for (int i = 0; i < numst; i++) {
    mt_Ss.push_back(mt_S);
    cholesky_decomp(stpdf[i].sigma,
                    mt_Ss[i],
                    diag);
  }

  std::vector<int> sample_per_gaussian(numst, 0);
  int samples_cnt = 0;
  for (int i = 0; i < numst - 1; i++) {
    sample_per_gaussian[i] = a00[i] * n;
    samples_cnt += sample_per_gaussian[i];
  }
  sample_per_gaussian[numst-1] = n - samples_cnt;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) {
      seq[i*dim+j] = gauss_std(gen);
      //      ++hist[std::round(seq[i*dim+j])];
    }
  }
  for (int i = 0, start_idx = 0; i < numst; i++) {  // i is the idx of state
    gauss_sample(seq, start_idx, sample_per_gaussian[i], i, mt_Ss[i], diag);
    start_idx += sample_per_gaussian[i];
  }
}
