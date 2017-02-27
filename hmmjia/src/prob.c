#include <math.h>
#include <stdio.h>
#include "hmm.h"
#define Pi 3.141592653589793
#define LOG_2_PI 1.83787706640935

int newgauss(GaussModel *md, int dim, int cls, int exist)
{
  int i,j,k,m,n;

  md->dim=dim;
  md->cls=cls;
  md->exist=exist;

  md->mean=(double *)calloc(dim,sizeof(double));
  matrix_2d_double(&(md->sigma),dim,dim);
  matrix_2d_double(&(md->sigma_inv), dim, dim);
  return 0;
}

int cpgauss(GaussModel *md1, GaussModel *md2)
{
  int i,j,k,m,n;

  md2->dim=md1->dim;
  md2->cls=md1->cls;
  md2->exist=md1->exist;
  md2->sigma_det=md1->sigma_det;

  for (i=0;i<md1->dim;i++)
    md2->mean[i]=md1->mean[i];
  
  for (i=0;i<md1->dim;i++){
    for (j=0;j<md1->dim;j++) {
      md2->sigma[i][j]=md1->sigma[i][j];
      md2->sigma_inv[i][j]=md1->sigma_inv[i][j];
    }
  }
  return 0;
}


void newhmm(HmmModel *md, int dim, int numst, int numcls, int *stcls)
{
  int i,j,k,m,n;
  int cls, exist=1;

  md->dim=dim;
  md->numst=numst;
  md->numcls=numcls;
  
  md->stcls=(int *)calloc(numst,sizeof(int));
  if (numcls==numst || stcls==NULL) {
    for (i=0;i<numst;i++) md->stcls[i]=i;
  }
  else {
    for (i=0;i<numst;i++) md->stcls[i]=stcls[i];
  }

  md->stpdf=(GaussModel **)calloc(numst, sizeof(GaussModel *));
  for (i=0; i<numst; i++) {
    md->stpdf[i]=(GaussModel *)calloc(1,sizeof(GaussModel));
    newgauss(md->stpdf[i], dim, md->stcls[i], exist);
  }
  
  matrix_2d_double(&md->a, numst,numst);
  md->a00=(double *)calloc(numst, sizeof(double));
}


void freehmm(HmmModel **md_pt)
{
  int i,j,k,m,n;
  int cls, exist=1,numst;
  HmmModel *md;

  md= *md_pt;
  numst=md->numst;

  for (i=0; i<numst; i++) {
    free(md->stpdf[i]->mean);
    free_matrix_2d_double(&(md->stpdf[i]->sigma), md->dim);
    free_matrix_2d_double(&(md->stpdf[i]->sigma_inv), md->dim);
    free(md->stpdf[i]);
  }
  free(md->stpdf);

  free(md->a00);
  free_matrix_2d_double(&md->a, numst);

  free(md->stcls);

  free(md);
  *md_pt=NULL;
}


void cphmm(HmmModel *md1, HmmModel *md2)
{
  int i,j,k,m,n;
  int cls, exist=1;
  int numst, numcls,dim;

  md2->dim=dim=md1->dim;
  md2->numst=numst=md1->numst;
  md2->numcls=numcls=md1->numcls;
  
  for (i=0;i<numst;i++) md2->stcls[i]=md1->stcls[i];

  for (i=0; i<numst; i++) {
    cpgauss(md1->stpdf[i], md2->stpdf[i]);
  }

  for (i=0;i<numst;i++) md2->a00[i]=md1->a00[i];
  for (i=0;i<numst;i++)
    for (j=0;j<numst;j++)
      md2->a[i][j]=md1->a[i][j];
}


double gauss_pdf_log(float *ft, GaussModel *gm)
{
  double res, tpdb, tpdb2, *db_array, *dif, *ptrdb1, *ptrdb2, *ptrdb3;
  int i,j,k,m,n;
  float *ptrft;

  if (!vector_double(&db_array, gm->dim))
    exit(1);

  if (!vector_double(&dif, gm->dim))
    exit(1);

  m=gm->dim;
  ptrdb1 = dif;
  ptrft = ft;
  ptrdb3 = gm->mean;
  for (i=0; i<m; i++) {
    *(ptrdb1++) = (double)(*ptrft)-(*ptrdb3);
    ptrft++;
    ptrdb3++;
  }

  ptrdb1 = db_array;
  for (i=0; i<m; i++) {
    *ptrdb1 = 0.0;
    ptrdb2 = gm->sigma_inv[i];
    ptrdb3 = dif;
    for (j=0; j<m; j++) {
      (*ptrdb1) += (*ptrdb2)*(*ptrdb3);
      ptrdb2++;
      ptrdb3++;
    }
    ptrdb1++;
  }

  tpdb2 = 0.0;
  ptrdb1 = db_array;
  ptrdb2 = dif;
  for (i=0; i<m; i++)
    {
      tpdb2 += (*ptrdb1)*(*ptrdb2);
      ptrdb1++;
      ptrdb2++;
    }

  tpdb = -((double)(gm->dim))/2.0*LOG_2_PI-0.5*log(gm->sigma_det);
  res = tpdb +(-0.5)*tpdb2;

  free(db_array);
  free(dif);

  return(res);
}

double gauss_pdf(float *ft, GaussModel *gm)
{
  return(exp(gauss_pdf_log(ft, gm)));
}

double mix_gauss_pdf_log(float *ft, GaussModel **gmlist, double *prior, 
			 int ncmp)
{
  double res, *h, v1,v2;
  int i,j,k,m,n;

  h=(double *)calloc(ncmp,sizeof(double));
  for (i=0;i<ncmp;i++)
    h[i]=gauss_pdf_log(ft, gmlist[i]);

  v1=h[0];
  for (i=1;i<ncmp;i++) if (h[i]>v1) v1=h[i];

  for (i=0,v2=0.0;i<ncmp;i++) {
    v2 += prior[i]*exp(h[i]-v1);
  }

  if (v2>0.0)  res=v1+log(v2); else res=-HUGE;

  free(h);
  return(res);
}

