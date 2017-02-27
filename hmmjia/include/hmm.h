#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "cluster.h"
#define EPSILON 1.0e-5
#define LAMBDA 0.999

typedef struct gaussmodel_struct
{
  int dim;
  int cls;
  int exist;
  double *mean;
  double **sigma;
  double **sigma_inv;
  double sigma_det;
} GaussModel;

typedef struct hmmmodel_struct
{
  int dim;
  int numst;   /* numst is the number of states */
  int numcls;
  int *stcls;
  GaussModel **stpdf;
  double **a;
  double *a00;  /* pmf of states at the boundary when there's no neighbor */
} HmmModel;


/*---------- functions ---------------*/
/*--------- estimate.c ---------------*/
/*------------------------------------*/

extern void forward(float *u, int ncols, double *thetalog, 
		    HmmModel *md, double *loglikehd);
extern void backward(float *u, int ncols, double *betalog, HmmModel *md);
extern void CompLm(float *u, int ncols, double *thetalog, double *betalog, 
		   double *Lm, HmmModel *md);
extern void CompHml(float *u, int ncols, double *thetalog, double *betalog, 
		    double *Hml, int west, HmmModel *md);
extern void viterbi(HmmModel *md, float *u, int len, int *optst, 
		    double *inita, double *lastmerit);
extern void viterbicls(HmmModel *md, float *u, int len, int *optst, double *inita,
		       double *lastmerit, int *bestnext);
extern void viterbi_mulseq(HmmModel *md, float **u, int nseq, int *len, 
			   int **st);
extern void updatepar_adder(float *u, int ncols, double *thetalog, 
			    double *betalog, double loglikehd, HmmModel *md, 
			    double *musum, double **mom2sum, double **asum, 
			    double *lsum);
extern void initialize(float **u, int nseq, int *len, int dim, HmmModel *md,
		       int ranflag);
extern double comploglike(HmmModel *md, float **u, int nseq, int *len, 
			  double *wt);
extern double classlikehd(HmmModel *md, float **u, int nseq, int *len, 
			  double **cprob, double *wt);
extern int baumwelch(float **u, int nseq, int *len, HmmModel *md, 
		      double *loglikehd, double *lhsumpt, double epsilon, 
		      double *wt);
extern void hmmfit(float **u, int nseq, int *len, int dim, HmmModel *md, 
		   int numst, int *stcls, double *loglikehd, double *lhsumpt, 
		   double epsilon, double *wt);

/*-------------------------------------*/
/*-------------- prob.c ---------------*/
/*-------------------------------------*/

extern int newgauss(GaussModel *md, int dim, int cls, int exist);
extern int cpgauss(GaussModel *md1, GaussModel *md2);
extern void newhmm(HmmModel *md, int dim, int numst, int numcls, int *stcls);
extern void freehmm(HmmModel **md_pt);
extern void cphmm(HmmModel *md1, HmmModel *md2);
extern double gauss_pdf_log(float *ft, GaussModel *gm);
extern double gauss_pdf(float *ft, GaussModel *gm);
extern double mix_gauss_pdf_log(float *ft, GaussModel **gmlist, double *prior, int ncmp);

/*-------------------------------------*/
/*------------- modelio.c -------------*/
/*-------------------------------------*/

extern unsigned char write_model(HmmModel *md, FILE *outfile);
extern unsigned char read_model(HmmModel *md, FILE *infile);
extern unsigned char print_model(HmmModel *md, FILE *outfile);


