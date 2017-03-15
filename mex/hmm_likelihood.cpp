//
//  hmm_likelihood_jia.c
//
//  Created by Yukun Chen on 11/20/15.
//  Copyright © 2015 Yukun Chen. All rights reserved.
//

#include "hmm.h"
#include "utils.h"
#include "mex.h"
#include <string.h>

#define MEX_ARGS int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs

// Do CHECK and throw a Mex error if check fails
inline void mxCHECK(bool expr, const char* msg) {
  if (!expr) {
    mexErrMsgTxt(msg);
  }
}
inline void mxERROR(const char* msg) { mexErrMsgTxt(msg); }


void mexPrint_mat_float(float* X, int dim, int n){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            mexPrintf("%f,\t", X[j+i*dim]);
        }
            mexPrintf(";\n");
    }
}

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

unsigned char mexPrint_model(HmmModel *md)
{
  int i,j,k,m,n;
  int dim,numst,numcls;

  dim=md->dim;
  numst=md->numst;
  numcls=md->numcls;
  
  mexPrintf("dim=%d\n", dim);
  mexPrintf("numst=%d\n", numst);
  mexPrintf("numcls=%d\n", numcls);

  mexPrintf("\nState class membership:\n");
  for (i=0; i<numst; i++)
    mexPrintf("%d ", md->stcls[i]);
  mexPrintf("\n");

  mexPrintf("\nTransition probability a00:\n");
  for (i=0; i<numst; i++)
    mexPrintf("%8e ", md->a00[i]);
  mexPrintf("\n");

  mexPrintf("Transition probability a:\n");
  for (i=0; i<numst; i++) {
    for (j=0; j<numst; j++)
      mexPrintf("%8e ", md->a[i][j]);
    mexPrintf("\n");
  }

  mexPrintf("\nThe Gaussian distributions of states:\n");
  for (i=0; i<numst; i++) {
    mexPrintf("\nState %d =============================\n", i);
    mexPrintf("exist=%d, dim=%d\n", md->stpdf[i]->exist, 
      md->stpdf[i]->dim);

    mexPrintf("Mean vector:\n");
    for (j=0; j<dim; j++)
      mexPrintf("%.5e ", md->stpdf[i]->mean[j]);
    mexPrintf("\n");

    mexPrintf("Sigma_det=%e\n",md->stpdf[i]->sigma_det);

    mexPrintf("Covariance matrix Sigma:\n");
 
    for (m=0; m<md->stpdf[i]->dim; m++) {
      for (n=0; n<md->stpdf[i]->dim; n++)
      mexPrintf("%.5e ", md->stpdf[i]->sigma[m][n]);
      mexPrintf("\n");
    }

    mexPrintf("Covariance matrix inverse Sigma_inv:\n");
 
    for (m=0; m<md->stpdf[i]->dim; m++) {
      for (n=0; n<md->stpdf[i]->dim; n++)
      mexPrintf("%.5e ", md->stpdf[i]->sigma_inv[m][n]);
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
  float *dat;
  float **u;
  double *wt=NULL;
  int nseq, numdat, onelen=0, *len, *stcls;
  HmmModel *md=NULL;
  int numst=2;
  double *loglikehd, lhsum;
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
  // size_t buflen;
  // int status;
  // if (!mxIsChar(prhs[0]) || (mxGetM(prhs[0]) != 1 ) )  {
  //     mexErrMsgIdAndTxt( "MATLAB:mxmalloc:invalidInput", 
  //               "Input infilename must be a string.");
  // }
  // /* Get number of characters in the input string.  Allocate enough
  //  memory to hold the converted string. */
  // buflen = mxGetN(prhs[0]) + 1;
  // infilename = mxMalloc(buflen);
  // /* Copy the string data into buf. */ 
  // status = mxGetString(prhs[0], infilename, (mwSize)buflen);   
  // mexPrintf("The input string is:  %s\n", infilename);

  // if (!mxIsChar(prhs[1]) || (mxGetM(prhs[1]) != 1 ) )  {
  //     mexErrMsgIdAndTxt( "MATLAB:mxmalloc:invalidInput", 
  //               "Input mdfilename must be a string.");
  // }
  // /* Get number of characters in the input string.  Allocate enough
  //  memory to hold the converted string. */
  // buflen = mxGetN(prhs[1]) + 1;
  // mdfilename = mxMalloc(buflen);
  // /* Copy the string data into buf. */ 
  // status = mxGetString(prhs[1], mdfilename, (mwSize)buflen);   
  // mexPrintf("The input string is:  %s\n", mdfilename);

  dim = mxGetScalar(prhs[1]);
  nseq = mxGetScalar(prhs[2]);
  onelen = mxGetScalar(prhs[3]);
  numst = mxGetScalar(prhs[4]);
  int verbose = mxGetScalar(prhs[5]);

  /*----------------------------------------------------------------*/
  /*----------------- Read in data ---------------------------------*/
  /*----------------------------------------------------------------*/

  // Assume the same length for all the sequences
  len = (int *)calloc(nseq, sizeof(int));
  if (onelen>0) {
    for (i=0;i<nseq;i++) len[i]=onelen;
  }
  else {
    for (i=0;i<nseq;i++) {     //read in len from stdin
      fscanf(stdin, "%d", len+i);
    }
  }

  for (i=0,numdat=0;i<nseq;i++) { numdat+=len[i];}
  double * dat_temp;
  dat_temp=(double *)mxGetData(prhs[0]);
  dat=(float *)calloc(numdat*dim,sizeof(float));
  //convert double * dat_temp to float * dat:
  for (i=0;i<numdat*dim;i++){
    dat[i]=dat_temp[i];
  }

  // dat=(float *)calloc(numdat*dim,sizeof(float));
  u=(float **)calloc(nseq,sizeof(float *));
  for (i=0,m=0;i<nseq;i++) {
    u[i]=dat+m*dim;
    m+=len[i];
  }

  int dimension=dim;
  int seq_len=onelen;
  if(verbose){
  mexPrintf("dimension=%d, n=%d\n",dimension,seq_len);
  mexPrintf("data:\n");
  mexPrint_mat_float(dat,dim,numdat);
  }

  /*----------------------------------------------------------------*/
  /*----------------- Estimate HMM  ---------------------------------*/
  /*----------------------------------------------------------------*/

  //fprintf(stderr, "numdat=%d, nseq=%d\n",numdat,nseq);

  // loglikehd=(double *)calloc(nseq,sizeof(double));
  

  // hmmfit(u, nseq, len, dim, md, numst, NULL, loglikehd, &lhsum,
	 // (double)epsilon, wt);


  /*-----------------------------------------------------------------*/
  /*----------------- Read in HMM   ---------------------------------*/
  /*-----------------------------------------------------------------*/
  md=(HmmModel *)calloc(1,sizeof(HmmModel));
  md->dim=dimension;
  md->numst=numst;
  md->numcls=numst;

  // prhs[6] = mxCreateDoubleMatrix(1, numst, mxREAL);
  double *pa00;
  pa00 = mxGetPr(prhs[6]);
  md->a00=(double *)calloc(numst,sizeof(double));
  for (i=0; i<md->numst; i++)
    md->a00[i] = *(pa00+i);
  if(verbose){
    mexPrintf("a00:\n");
    mexPrint_mat_double(md->a00,numst,1);
  }
  // prhs[1] = mxCreateDoubleMatrix(numst, numst, mxREAL);
  double *pa;
  pa = mxGetPr(prhs[7]);
  md->a=(double **)calloc(numst,sizeof(double *));
  for(i=0; i<md->numst; i++){
    md->a[i]=(double *)calloc(numst,sizeof(double));
    for(j=0; j<md->numst; j++){
      md->a[i][j]=*(pa+i*md->numst+j);
    }
  }
  if(verbose){
    mexPrintf("a:\n");
    for(i=0; i<md->numst; i++){
      mexPrint_mat_double(md->a[i],numst,1);
    }
  }
  // const int NUMBER_OF_FIELDS=2;

  /*----------- Assign State class membership ---------------------*/
  md->stcls=(int *)calloc(numst,sizeof(int));
  for (i=0; i<numst; i++)
    md->stcls[i]=i;

  md->stpdf=(GaussModel **)calloc(numst, sizeof(GaussModel *));
  for (i=0; i<numst; i++)
    md->stpdf[i]=(GaussModel *)calloc(1, sizeof(GaussModel));

  mwSize sizebuf;
  for (i=0; i<numst; i++) {
      md->stpdf[i]->exist=1;
      md->stpdf[i]->dim=dimension;
      /*----------- read mean vector ---------------------*/
      md->stpdf[i]->mean=(double *)calloc(md->stpdf[i]->dim,sizeof(double));
      double * mean_ptr = mxGetPr(prhs[8]);
      for (j=0; j<md->stpdf[i]->dim; j++){
          md->stpdf[i]->mean[j]=*(mean_ptr+j+i*dimension);
      }
      if(verbose){
        mexPrintf("i=%d \t ",i);
        for (j=0; j<dimension; j++){
          mexPrintf("%f (%f)\t",md->stpdf[i]->mean[j],*(mean_ptr+j+i*dimension));
        }
        mexPrintf("\n");
      }
      /*----------- read sigma       ---------------------*/
      md->stpdf[i]->sigma_det=1;
      md->stpdf[i]->sigma=(double **)calloc(md->stpdf[i]->dim,sizeof(double *));
      double *sigma_ptr = mxGetPr(prhs[9]);
      for (m=0; m<md->stpdf[i]->dim; m++) {
        md->stpdf[i]->sigma[m]=(double *)calloc(md->stpdf[i]->dim, 
            sizeof(double));
        for (n=0; n<md->stpdf[i]->dim; n++){
          md->stpdf[i]->sigma[m][n] = *(sigma_ptr+i*(md->stpdf[i]->dim*md->stpdf[i]->dim)+m*md->stpdf[i]->dim+n);
          }
      }
      if(verbose){
        mexPrintf("i=%d \t ",i);
        for(m=0; m<md->stpdf[i]->dim;m++){
          mexPrint_mat_double(md->stpdf[i]->sigma[m],md->stpdf[i]->dim,1);
        }
        mexPrintf("\n");
      }
      /*----------- read sigma_inv       ---------------------*/
      md->stpdf[i]->sigma_inv=(double **)calloc(md->stpdf[i]->dim,sizeof(double *));
      double *sigma_inv_ptr = mxGetPr(prhs[10]);
      for (m=0; m<md->stpdf[i]->dim; m++) {
        md->stpdf[i]->sigma_inv[m]=(double *)calloc(md->stpdf[i]->dim, 
            sizeof(double));
        for (n=0; n<md->stpdf[i]->dim; n++){
          md->stpdf[i]->sigma_inv[m][n] = *(sigma_inv_ptr+i*(md->stpdf[i]->dim*md->stpdf[i]->dim)+m*md->stpdf[i]->dim+n);
          }
      }
      if(verbose){
        mexPrintf("\n i=%d \t ",i);
        for(m=0; m<md->stpdf[i]->dim;m++){
          mexPrint_mat_double(md->stpdf[i]->sigma_inv[m],md->stpdf[i]->dim,1);
        }
        mexPrintf("\n");
      }
  }

  if(verbose){
    mexPrint_model(md);
  }
  /*----------------------------------------------------------------*/
  /*------------ Compute the likelihood of the sequence  -----------*/
  /*----------------------------------------------------------------*/
  // double loglikelihood,oneseq_likelihood;
  // double *thetalog;
  // thetalog=(double *)calloc(seq_len*md->numst, sizeof(double));
  // for (i=0, loglikelihood=0.0; i<nseq; i++) {
  //     forward(u[i], len[i], thetalog,  md, &oneseq_likelihood);
  //     mexPrintf("loglikelihood for %d seq = %f\n", i, oneseq_likelihood);
  //     loglikelihood+=oneseq_likelihood;
  // }

  double loglikelihood,oneseq_likelihood;
  double *thetalog;
  thetalog=(double *)calloc(seq_len*md->numst, sizeof(double));
  for (i=0, loglikelihood=0.0; i<nseq; i++) {
      forward(u[i], len[i], thetalog,  md, &oneseq_likelihood);
      // printf("loglikelihood for %d seq = %f\n", i, oneseq_likelihood);
      loglikelihood+=oneseq_likelihood;
  }
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *b;
  b = mxGetPr(plhs[0]);
  *b=loglikelihood; // /nseq/onelen

  free(dat);
  free(thetalog);
  freehmm(&md);

}
