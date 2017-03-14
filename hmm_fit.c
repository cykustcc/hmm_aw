//
//  hmm_fit_jia.c
//
//  Created by Yukun Chen on 10/20/15.
//  Copyright © 2015 Yukun Chen. All rights reserved.
//


#include "hmm.h"
#include "utils.h"
#include "mex.h"

void mexPrint_mat(float* X, int dim, int n){
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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
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

  if (nrhs != 6){
    mexErrMsgTxt("Must have 6 input arguments");
  }
  if (nlhs != 3){
    mexErrMsgTxt("Must have 3 output arguments");
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
  dat_temp=mxGetData(prhs[0]);
  dat=(float *)calloc(numdat*dim,sizeof(float));
  //convert double * dat_temp to float * dat:
  for (i=0;i<numdat*dim;i++){
    dat[i]=dat_temp[i];
  }

  // dat=(float *)calloc(numdat*dim,sizeof(float));
  u=(double **)calloc(nseq,sizeof(double *));
  for (i=0,m=0;i<nseq;i++) {
    u[i]=dat+m*dim;
    m+=len[i];
  }

  int dimension=dim;
  int seq_len=onelen;
  if(verbose){
  mexPrintf("dimension=%d, n=%d\n",dimension,seq_len);
  mexPrint_mat(dat,dim,numdat);
  }

  /*----------------------------------------------------------------*/
  /*----------------- Estimate HMM  ---------------------------------*/
  /*----------------------------------------------------------------*/

  //fprintf(stderr, "numdat=%d, nseq=%d\n",numdat,nseq);

  loglikehd=(double *)calloc(nseq,sizeof(double));
  md=(HmmModel *)calloc(1,sizeof(HmmModel));

  hmmfit(u, nseq, len, dim, md, numst, NULL, loglikehd, &lhsum,
	 (double)epsilon, wt);

  //Output loglikehd from hmmfit() is not written out

  // Binary file for the output model
  // write_model(md, mdfile);

  //Ascii file for the model
  if(verbose){
    mexPrint_model(md);
  }
  /*----------------------------------------------------------------*/
  /*----------------- Return HMM   ---------------------------------*/
  /*----------------------------------------------------------------*/
  plhs[0] = mxCreateDoubleMatrix(1, md->numst, mxREAL);
  double *pa00;
  pa00 = mxGetPr(plhs[0]);
  for (i=0; i<md->numst; i++)
    pa00[i] = md->a00[i];

  plhs[1] = mxCreateDoubleMatrix(md->numst, md->numst, mxREAL);
  double *pa;
  int count=0;
  pa = mxGetPr(plhs[1]);
  for(i=0; i<md->numst; i++){
    for(j=0; j<md->numst; j++){
      *(pa+count)=md->a[j][i];
      count++;
    }
  }
  const int NUMBER_OF_FIELDS=2;
  const int NUMBER_OF_STRUCTS=md->numst;
  mwSize dims[2] = {1, NUMBER_OF_STRUCTS };
  const char *field_names[] = {"means", "sigma"};
  plhs[2] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);
  int means_field = mxGetFieldNumber(plhs[2],"means");
  int sigma_field = mxGetFieldNumber(plhs[2],"sigma");
    /* Populate the name and phone fields of the phonebook structure. */ 
  for (i=0; i<NUMBER_OF_STRUCTS; i++) {
      /*----------- return mean vector ---------------------*/
      mxArray *field_mean_vector;
      field_mean_vector = mxCreateDoubleMatrix(1,dimension,mxREAL);
      double *p_field_mean_vector = mxGetPr(field_mean_vector);
      for (j=0; j<dimension; j++)
        p_field_mean_vector[j]=md->stpdf[i]->mean[j];
      mxSetFieldByNumber(plhs[2],i,means_field,field_mean_vector);
      /*----------- return sigma       ---------------------*/
      mxArray *field_sigma;
      field_sigma = mxCreateDoubleMatrix(dimension,dimension,mxREAL);
      double *p_field_sigma = mxGetPr(field_sigma);
      int count=0;
      for (m=0; m<md->stpdf[i]->dim; m++) {
        for (n=0; n<md->stpdf[i]->dim; n++){
          p_field_sigma[count] = md->stpdf[i]->sigma[n][m];
          count++;
          }
      }
      mxSetFieldByNumber(plhs[2],i,sigma_field,field_sigma);
  }
}
