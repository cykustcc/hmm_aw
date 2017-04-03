//
//  hmm_likelihood.c
//
//  Created by Yukun Chen on 11/20/15.
//  Copyright © 2015 Yukun Chen. All rights reserved.
//

#include "hmm.h"
#include "utils.h"
#include "mex.h"
#include <string.h>
#include <vector>

#define MEX_ARGS int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs

// Do CHECK and throw a Mex error if check fails
inline void mxCHECK(bool expr, const char* msg) {
  if (!expr) {
    mexErrMsgTxt(msg);
  }
}
inline void mxERROR(const char* msg) { mexErrMsgTxt(msg); }


void mexPrint_mat_double(double* X, int dim, int n){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            mexPrintf("%f,\t", X[j+i*dim]);
        }
            mexPrintf(";\n");
    }
}

template<typename DType>
void mex_print_vector(std::vector<DType>& vt,
                      /* # of rows*/ int m,
                      /* # of cols*/ int n){
  int size = vt.size();
  assert(m*n <= size);
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      mexPrintf("%f,\t", vt[i*n + j]);
    }
    mexPrintf(";\n");
  }
}

template<typename DType>
void mex_print_matrix(std::vector<std::vector<DType>> &mt){
  int rows = mt.size();
  int cols = rows > 0 ? mt[0].size(): 0;

  for (int i=0; i<rows; i++) {
    for (int j=0; j<cols; j++) {
      mexPrintf("%f,\t", mt[i][j]);
    }
    mexPrintf(";\n");
  }
}

unsigned char mexPrint_model(HmmModel &md)
{
  int i,j,k,m,n;
  int dim,numst,numcls;

  dim=md.dim;
  numst=md.numst;
  numcls=md.numcls;

  mexPrintf("dim=%d\n", dim);
  mexPrintf("numst=%d\n", numst);
  mexPrintf("numcls=%d\n", numcls);

  mexPrintf("\nState class membership:\n");
  for (i=0; i<numst; i++)
    mexPrintf("%d ", md.stcls[i]);
  mexPrintf("\n");

  mexPrintf("\nTransition probability a00:\n");
  for (i=0; i<numst; i++)
    mexPrintf("%8e ", md.a00[i]);
  mexPrintf("\n");

  mexPrintf("Transition probability a:\n");
  for (i=0; i<numst; i++) {
    for (j=0; j<numst; j++)
      mexPrintf("%8e ", md.a[i][j]);
    mexPrintf("\n");
  }

  mexPrintf("\nThe Gaussian distributions of states:\n");
  for (i=0; i<numst; i++) {
    mexPrintf("\nState %d =============================\n", i);
    mexPrintf("exist=%d, dim=%d\n", md.stpdf[i].exist,
      md.stpdf[i].dim);

    mexPrintf("Mean vector:\n");
    for (j=0; j<dim; j++)
      mexPrintf("%.5e ", md.stpdf[i].mean[j]);
    mexPrintf("\n");

    mexPrintf("Sigma_det=%e\n",md.stpdf[i].sigma_det);

    mexPrintf("Covariance matrix Sigma:\n");

    for (m=0; m<md.stpdf[i].dim; m++) {
      for (n=0; n<md.stpdf[i].dim; n++)
      mexPrintf("%.5e ", md.stpdf[i].sigma[m][n]);
      mexPrintf("\n");
    }

    mexPrintf("Covariance matrix inverse Sigma_inv:\n");

    for (m=0; m<md.stpdf[i].dim; m++) {
      for (n=0; n<md.stpdf[i].dim; n++)
      mexPrintf("%.5e ", md.stpdf[i].sigma_inv[m][n]);
      mexPrintf("\n");
    }
  }

  return(1);
}
void mexFunction(MEX_ARGS)
{

  char *infilename;
  char *mdfilename;
  FILE *infile, *mdfile;
  int i,j,k,m,n;
  int dim=2;
  int nseq, numdat, onelen=0;
  int numst=2;
  double lhsum;
  float epsilon=EPSILON;
  float tp1, tp2;

  // mexPrintf("nrhs = %d\n", nrhs);

  if (nrhs != 11){
    mexErrMsgTxt("Must have 11 input arguments");
  }
  if (nlhs != 1){
    mexErrMsgTxt("Must have 1 output arguments");
  }
  /*----------------------------------------------------------------*/
  /*---------------- Read in parameters from command line-----------*/
  /*----------------------------------------------------------------*/

  dim = mxGetScalar(prhs[1]);
  nseq = mxGetScalar(prhs[2]);
  onelen = mxGetScalar(prhs[3]);
  numst = mxGetScalar(prhs[4]);
  int verbose = mxGetScalar(prhs[5]);

  /*----------------------------------------------------------------*/
  /*----------------- Read in data ---------------------------------*/
  /*----------------------------------------------------------------*/

  // Assume the same length for all the sequences
  std::vector<int> len(nseq, 0);

  if (onelen>0) {
    for (i=0;i<nseq;i++) len[i]=onelen;
  }
  else {
    for (i=0;i<nseq;i++) {     //read in len from stdin
      fscanf(stdin, "%d", &len[i]);
    }
  }
  int maxlen = *std::max_element(len.begin(), len.end());

  for (i=0,numdat=0;i<nseq;i++) { numdat+=len[i];}

  double* dat_temp;
  dat_temp= (double *)mxGetData(prhs[0]);
  std::vector<std::vector<float>> u(nseq, std::vector<float>());

  for (i=0;i<nseq;i++) {
    for (int j = 0; j < dim*len[i]; j++)
    {
      u[i].push_back(dat_temp[j]);
    }
  }

  int dimension=dim;
  int seq_len=onelen;
  if(verbose){
    mexPrintf("dimension=%d, n=%d\n",dimension,seq_len);
    for (i=0;i<nseq;i++) {
      mex_print_vector(u[i],len[i],dim);
    }
  }

  /*-----------------------------------------------------------------*/
  /*----------------- Read in HMM   ---------------------------------*/
  /*-----------------------------------------------------------------*/
  std::vector<int> stcls;

  HmmModel md(dimension, numst, numst, stcls);

  // prhs[6] = mxCreateDoubleMatrix(1, numst, mxREAL);
  double *pa00;
  pa00 = mxGetPr(prhs[6]);
  for (i=0; i<md.numst; i++)
    md.a00[i] = *(pa00+i);
  if(verbose){
    mexPrintf("a00:\n");
    // mexPrint_mat_double(md.a00,numst,1);
    mex_print_vector(md.a00,1,numst);
  }
  // prhs[1] = mxCreateDoubleMatrix(numst, numst, mxREAL);
  double *pa;
  pa = mxGetPr(prhs[7]);
  for(i=0; i<md.numst; i++){
    for(j=0; j<md.numst; j++){
      md.a[i][j]=*(pa+i*md.numst+j);
    }
  }
  if(verbose){
    mexPrintf("a:\n");
    for(i=0; i<md.numst; i++){
      // mexPrint_mat_double(md.a[i],numst,1);
      mex_print_vector(md.a[i],1,numst);
    }
  }
  // const int NUMBER_OF_FIELDS=2;

  /*----------- Assign State class membership ---------------------*/

  mwSize sizebuf;
  for (i=0; i<numst; i++) {
      md.stpdf[i].exist=1;
      md.stpdf[i].dim=dimension;
      /*----------- read mean vector ---------------------*/
      double * mean_ptr = mxGetPr(prhs[8]);
      for (j=0; j<md.stpdf[i].dim; j++){
          md.stpdf[i].mean[j]=*(mean_ptr+j+i*dimension);
      }
      if(verbose){
        mexPrintf("i=%d \t ",i);
        for (j=0; j<dimension; j++){
          mexPrintf("%f (%f)\t",md.stpdf[i].mean[j],*(mean_ptr+j+i*dimension));
        }
        mexPrintf("\n");
      }
      /*----------- read sigma       ---------------------*/
      md.stpdf[i].sigma_det=1;
      double *sigma_ptr = mxGetPr(prhs[9]);
      for (m=0; m<md.stpdf[i].dim; m++) {
        for (n=0; n<md.stpdf[i].dim; n++){
          md.stpdf[i].sigma[m][n] = *(sigma_ptr+i*(md.stpdf[i].dim*md.stpdf[i].dim)+m*md.stpdf[i].dim+n);
          }
      }
      if(verbose){
        mexPrintf("i=%d \t ",i);
        mex_print_matrix(md.stpdf[i].sigma);
        // for(m=0; m<md.stpdf[i].dim;m++){
        //   mexPrint_mat_double(md.stpdf[i].sigma[m],md.stpdf[i].dim,1);
        // }
        mexPrintf("\n");
      }
      /*----------- read sigma_inv       ---------------------*/
      double *sigma_inv_ptr = mxGetPr(prhs[10]);
      for (m=0; m<md.stpdf[i].dim; m++) {
        for (n=0; n<md.stpdf[i].dim; n++){
          md.stpdf[i].sigma_inv[m][n] = *(sigma_inv_ptr+i*(md.stpdf[i].dim*md.stpdf[i].dim)+m*md.stpdf[i].dim+n);
          }
      }
      int tmp=mat_det_inv(md.stpdf[i].sigma, md.stpdf[i].sigma_inv,
                              md.stpdf[i].sigma_det,md.stpdf[i].dim);
      if(verbose){
        mexPrintf("\n i=%d \t ",i);
        mex_print_matrix(md.stpdf[i].sigma_inv);
        // for(m=0; m<md.stpdf[i].dim;m++){
        //   mexPrint_mat_double(md.stpdf[i].sigma_inv[m],md.stpdf[i].dim,1);
        // }
        mexPrintf("\n");
      }
  }

  if(verbose){
    mexPrint_model(md);
  }
  /*----------------------------------------------------------------*/
  /*------------ Compute the likelihood of the sequence  -----------*/
  /*----------------------------------------------------------------*/

  double loglikelihood,oneseq_likelihood;
  std::vector<double> thetalog(seq_len*md.numst, 0.0);
  std::vector<double> wt;
  // loglikelihood = md.comploglike(u, 1, len, wt);
  for (i=0, loglikelihood=0.0; i<nseq; i++) {

      md.forward(u[i], len[i], thetalog, oneseq_likelihood);
      // printf("loglikelihood for %d seq = %f\n", i, oneseq_likelihood);
      loglikelihood += oneseq_likelihood;
      if (verbose)
        mex_print_vector(thetalog, len[i], numst);
  }

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *b;
  b = mxGetPr(plhs[0]);
  *b=loglikelihood; // /nseq/onelen

}
