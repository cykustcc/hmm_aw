#ifndef __hmm_aw_hmm_h
#define __hmm_aw_hmm_h

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "cluster.h"
#include <vector>
#include <cassert>
#include <algorithm>

#define EPSILON 1.0e-5
#define LAMBDA 0.999

class GaussModel{
public:
  int dim;
  int cls;
  int exist;
  std::vector<double> mean;
  std::vector<std::vector<double>> sigma;
  std::vector<std::vector<double>> sigma_inv;
  double sigma_det;
  GaussModel(){
    GaussModel(1,1,1);
  }
  GaussModel(int _dim, int _cls, int _exist): dim(_dim), cls(_cls), exist(_exist){
    mean.resize(dim);
    sigma.resize(dim);
    for(int i = 0; i < dim; i++)
      sigma[i].resize(dim);
    sigma_inv.resize(dim);
    for(int i = 0; i < dim; i++)
      sigma_inv[i].resize(dim);
  }
  GaussModel(const GaussModel &g): GaussModel(g.dim, g.cls, g.exist){
    for (int i=0; i<g.dim; i++) {
      mean[i] = g.mean[i];
    }
    sigma_det=g.sigma_det;
    for (int i=0; i<g.dim; i++) {
      for (int j=0; j<g.dim; j++) {
        sigma[i][j] = g.sigma[i][j];
        sigma_inv[i][j] = g.sigma_inv[i][j];
      }
    }
  }
  GaussModel & operator=(const GaussModel &g){
    dim = g.dim;
    cls = g.cls;
    exist = g.exist;
    mean.resize(dim);
    for (int i=0; i<g.dim; i++) {
      mean[i] = g.mean[i];
    }
    sigma_det=g.sigma_det;
    sigma.resize(dim);
    sigma_inv.resize(dim);
    for (int i=0; i<g.dim; i++) {
      sigma[i].resize(dim);
      sigma_inv[i].resize(dim);
    }
    for (int i=0; i<g.dim; i++) {
      for (int j=0; j<g.dim; j++) {
        sigma[i][j] = g.sigma[i][j];
        sigma_inv[i][j] = g.sigma_inv[i][j];
      }
    }
    return *this;
  }
  ~GaussModel() = default;
  /*-------------------------------------*/
  /*-------------- prob.c ---------------*/
  /*-------------------------------------*/
  
  double gauss_pdf_log(std::vector<float> &ft,
                       int baseidx) const;
  
  double gauss_pdf(std::vector<float> &ft,
                   int baseidx) const;
};

extern double mix_gauss_pdf_log(std::vector<float> &ft,
                                std::vector<GaussModel> &gmlist,
                                std::vector<double> &prior,
                                int ncmp,
                                int baseidx);

class HmmModel{
public:
  int dim;
  int numst; /* numst is the number of states */
  int numcls;
  std::vector<int> stcls;
  std::vector<GaussModel> stpdf;
  std::vector<std::vector<double>> a;
  std::vector<double> a00; /* pmf of states at the boundary when there's no neighbor */
  HmmModel(){
    std::vector<int> stcls_tmp;
    HmmModel(0,0,0,stcls_tmp);
  }
  HmmModel(int _dim, int _numst, int _numcls, std::vector<int>& _stcls)
    : dim(_dim), numst(_numst), numcls(_numcls){
    assert(_stcls.size() == 0 || _stcls.size() == numst);
    stcls.resize(numst);
    if (numcls == numst || stcls.size() == 0) {
      for (int i=0; i<numst; i++) stcls[i]=i;
    }else{
      for (int i=0; i<numst; i++) stcls[i] = _stcls[i];
    }
//    stpdf.resize(numst);
    for (int i=0; i<numst; i++) {
      GaussModel tmp(dim, stcls[i], 1);
      stpdf.push_back(tmp);
    }
    a.resize(numst);
    for (int i=0; i<numst; i++) {
      a[i].resize(numst);
    }
    a00.resize(numst);
  }
  HmmModel(HmmModel &h): HmmModel(h.dim, h.numst, h.numcls, h.stcls){
    for (int i=0; i<numst; i++) {
      stcls[i] = h.stcls[i];
      stpdf[i] = h.stpdf[i];
      a00[i] = h.a00[i];
      for (int j=0; j<numst; j++) {
        a[i][j] = h.a[i][j];
      }
    }
  }
  
  HmmModel & operator=(HmmModel &h){
    dim = h.dim;
    numst = h.numst;
    numcls = h.numcls;
    stcls.resize(numst);
    a.resize(numst);
    for (int i=0; i<numst; i++) {
      a[i].resize(numst);
    }
    a00.resize(numst);
    stpdf.resize(numst);
    for (int i=0; i<numst; i++) {
      stcls[i] = h.stcls[i];
      stpdf[i] = h.stpdf[i];
      a00[i] = h.a00[i];
      for (int j=0; j<numst; j++) {
        a[i][j] = h.a[i][j];
      }
    }
    return *this;
  }
  ~HmmModel() = default;
  void resize(int _dim, int _numst, int _numcls, std::vector<int> _stcls){
    dim = _dim;
    numst = _numst;
    numcls = _numcls;
    assert(_stcls.size() == 0 || _stcls.size() == numst);
    stcls.resize(numst);
    if (numcls == numst || stcls.size() == 0) {
      for (int i=0; i<numst; i++) stcls[i]=i;
    }else{
      for (int i=0; i<numst; i++) stcls[i] = _stcls[i];
    }
    //    stpdf.resize(numst);
    for (int i=0; i<numst; i++) {
      GaussModel tmp(dim, stcls[i], 1);
      stpdf.push_back(tmp);
    }
    a.resize(numst);
    for (int i=0; i<numst; i++) {
      a[i].resize(numst);
    }
    a00.resize(numst);
  }
  /*--------- model io         ----------*/
  void read_model(std::string filename);
  void write_model(std::string filename) const;
  void print_model(std::string filename) const;
  /*--------- model estimation ----------*/
  /*--------- estimate.c ----------------*/
  void forward(std::vector<float> &u,
               int ncols,
               std::vector<double> &thetalog,
               double &loglikehd);
  void backward(std::vector<float> &u,
               int ncols,
               std::vector<double> &betalog);
  void CompLm(std::vector<float> &u,
              int ncols,
              std::vector<double> &thetalog,
              std::vector<double> &betalog,
              std::vector<double> &Lm);
  void CompHml(std::vector<float> &u,
               int ncols,
               std::vector<double> &thetalog,
               std::vector<double> &betalog,
               std::vector<double> &Hml,
               int west);
  void viterbi(std::vector<float> &u,
               int len, std::vector<int> &optst,
               std::vector<double> & inita,
               std::vector<double> &lastmerit);
  
  void formmix(std::vector<double> &inita,
               std::vector<std::vector<double>> &tm,
               std::vector<double> &astart,
               std::vector<std::vector<GaussModel>> &pdflist,
               std::vector<std::vector<double>> &prior,
               std::vector<int> &nstpercls);
  
  void viterbicls(std::vector<float> &u,
                  int len,
                  std::vector<int> &optst,
                  std::vector<double> &inita,
                  std::vector<double> &lastmerit,
                  int &bestnext);
  
  void viterbi_mulseq(std::vector<std::vector<float>> &u,
                      int nseq,
                      std::vector<int> &len,
                      std::vector<std::vector<int>> &st);
  
  void updatepar_adder(std::vector<float> &u,
                       int ncols,
                       std::vector<double> &thetalog,
                       std::vector<double> &betalog,
                       double loglikehd,
                       std::vector<double> &musum,
                       std::vector<std::vector<double>> &mom2sum,
                       std::vector<std::vector<double>> &asum,
                       std::vector<double> &lsum);
  
  void initialize(std::vector<std::vector<float>> &u,
                  int nseq,
                  std::vector<int> &len,
                  int dim,
                  int ranflag);
  
  double comploglike(std::vector<std::vector<float>> &u,
                     int nseq,
                     std::vector<int> &len,
                     std::vector<double>&wt);
  
  double classlikehd(std::vector<std::vector<float>> &u,
                     int nseq,
                     std::vector<int> &len,
                     std::vector<std::vector<double>> &cprob,
                     std::vector<double>& wt);
  
  int baumwelch(std::vector<std::vector<float>> &u,
                int nseq,
                std::vector<int> &len,
                std::vector<double> &loglikehd,
                double &lhsumpt,
                double epsilon,
                std::vector<double> &wt,
                bool forcediag);
  
  void hmmfit(std::vector<std::vector<float>> &u,
              int nseq,
              std::vector<int> &len,
              std::vector<double> &loglikehd,
              double &lhsumpt,
              double epsilon,
              std::vector<double> &wt,
              bool forcediag);
  
  double dist_KL(HmmModel &hmm2, int sample_size, bool diag = false);
  double dist_MAW(HmmModel &hmm2);
  double dist_IAW(HmmModel &hmm2);
  
  void gen_seq(std::vector<float> &seq,
               int n,
               bool diag);
  
  void gauss_sample(std::vector<float> &seq,
                    int idx,
                    int state,
                    std::vector<std::vector<double>> &mt_S,
                    bool diag);
  
};







#endif


