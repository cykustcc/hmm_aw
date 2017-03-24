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
  ~GaussModel(){
//    mean.clear();
//    sigma.clear();
//    sigma_inv.clear();
  }
};

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
  HmmModel(int _dim, int _numst, int _numcls, std::vector<int>& _stcls): dim(_dim), numst(_numst), numcls(_numcls){
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
  ~HmmModel(){
//    stcls.clear();
//    stpdf.clear();
//    a.clear();
//    a00.clear();
  }
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
  void read_model(std::string filename);
  void write_model(std::string filename);
  void print_model(std::string filename);
};


/*---------- functions ---------------*/
/*--------- estimate.c ---------------*/
/*------------------------------------*/

extern void forward(std::vector<float> &u,
                    int ncols,
                    std::vector<double> &thetalog,
                    HmmModel &md,
                    double &loglikehd);

extern void backward(std::vector<float> &u,
                     int ncols,
                     std::vector<double> &betalog,
                     HmmModel &md);

extern void CompLm(std::vector<float> &u,
                   int ncols,
                   std::vector<double> &thetalog,
                   std::vector<double> &betalog,
                   std::vector<double> &Lm,
                   HmmModel &md);

extern void CompHml(std::vector<float> &u,
                    int ncols,
                    std::vector<double> &thetalog,
                    std::vector<double> &betalog,
                    std::vector<double> &Hml,
                    int west,
                    HmmModel &md);

extern void viterbi(HmmModel &md,
                    std::vector<float> &u,
                    int len, std::vector<int> &optst,
                    std::vector<double> & inita,
                    std::vector<double> &lastmerit);

extern void viterbicls(HmmModel &md,
                       std::vector<float> &u,
                       int len,
                       std::vector<int> &optst,
                       std::vector<double> &inita,
                       std::vector<double> &lastmerit,
                       int &bestnext);

extern void viterbi_mulseq(HmmModel &md,
                           std::vector<std::vector<float>> &u,
                           int nseq,
                           std::vector<int> &len,
                           std::vector<std::vector<int>> &st);

extern void updatepar_adder(std::vector<float> &u,
                            int ncols,
                            std::vector<double> &thetalog,
                            std::vector<double> &betalog,
                            double loglikehd,
                            HmmModel &md,
                            std::vector<double> &musum,
                            std::vector<std::vector<double>> &mom2sum,
                            std::vector<std::vector<double>> &asum,
                            std::vector<double> &lsum);

extern void initialize(std::vector<std::vector<float>> &u,
                       int nseq,
                       std::vector<int> &len,
                       int dim,
                       HmmModel &md,
                       int ranflag);

extern double comploglike(HmmModel &md,
                          std::vector<std::vector<float>> &u,
                          int nseq,
                          std::vector<int> &len,
                          std::vector<double>&wt);

extern double classlikehd(HmmModel &md,
                          std::vector<std::vector<float>> &u,
                          int nseq,
                          std::vector<int> &len,
                          std::vector<std::vector<double>> &cprob,
                          std::vector<double>& wt);

extern int baumwelch(std::vector<std::vector<float>> &u,
                     int nseq,
                     std::vector<int> &len,
                     HmmModel &md,
                     std::vector<double> &loglikehd,
                     double &lhsumpt,
                     double epsilon,
                     std::vector<double> &wt,
                     bool forcediag);

extern void hmmfit(HmmModel& md,
                   std::vector<std::vector<float>> &u,
                   int nseq,
                   std::vector<int> &len,
                   std::vector<double> &loglikehd,
                   double &lhsumpt,
                   double epsilon,
                   std::vector<double> &wt,
                   bool forcediag);

/*-------------------------------------*/
/*-------------- prob.c ---------------*/
/*-------------------------------------*/

extern double gauss_pdf_log(std::vector<float> &ft,
                            GaussModel &gm,
                            int baseidx);

extern double gauss_pdf(std::vector<float> &ft,
                        GaussModel &gm,
                        int baseidx);

extern double mix_gauss_pdf_log(std::vector<float> &ft,
                                std::vector<GaussModel> &gmlist,
                                std::vector<double> &prior,
                                int ncmp,
                                int baseidx);

/*-------------------------------------*/
/*------------- modelio.c -------------*/
/*-------------------------------------*/

//extern unsigned char write_model(  HmmModel &md, FILE *outfile);
//extern unsigned char read_model(  HmmModel &md, FILE *infile);
//extern unsigned char print_model(  HmmModel &md, FILE *outfile);
//int hmm_read(  HmmModel &md, const char* filename);
//int hmm_write(  HmmModel &md, const char* filename);

#endif


