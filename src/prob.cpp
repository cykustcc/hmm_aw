#include <math.h>
#include <stdio.h>
#include "hmm.h"
#define Pi 3.141592653589793
#define LOG_2_PI 1.83787706640935

double GaussModel::gauss_pdf_log(std::vector<float> &ft,
                                 int baseidx) const {
  double res, tpdb, tpdb2;
  std::vector<double> db_array(dim, 0.0);
  std::vector<double> dif(dim, 0.0);

  int m = dim;
  for (int i = 0; i < m; i++) {
    dif[i] = ft[baseidx + i] - mean[i];
  }

//  ptrdb1 = db_array;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      db_array[i] += sigma_inv[i][j] * dif[j];
    }
  }

  tpdb2 = 0.0;
  for (int i = 0; i < m; i++) {
    tpdb2 += db_array[i] * dif[i];
  }

  tpdb = -((double)(dim)) / 2.0 * LOG_2_PI - 0.5 * log(sigma_det);
  res = tpdb + (-0.5) * tpdb2;

  return res;
}

double GaussModel::gauss_pdf(std::vector<float> &ft,
                             int baseidx) const {
  return(exp(gauss_pdf_log(ft, baseidx)));
}

double mix_gauss_pdf_log(std::vector<float> &ft,
                         std::vector<GaussModel> &gmlist,
                         std::vector<double> &prior,
                         int ncmp,
                         int baseidx) {
  double res, v1, v2;

  std::vector<double> h(ncmp, 0.0);
  for (int i = 0; i < ncmp; i++)
    h[i] = gmlist[i].gauss_pdf_log(ft, baseidx);

  v1 = h[0];
  for (int i = 1; i < ncmp; i++)
    if (h[i] > v1)
      v1 = h[i];

  v2 = 0.0;
  for (int i = 0; i < ncmp; i++) {
    v2 += prior[i] * exp(h[i] - v1);
  }

  if (v2 > 0.0)
    res = v1 + log(v2);
  else
    res = -HUGE;

  return(res);
}

