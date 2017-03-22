#include <math.h>
#include <stdio.h>
#include "hmm.h"
#define Pi 3.141592653589793
#define LOG_2_PI 1.83787706640935

//int newgauss(GaussModel *md, int dim, int cls, int exist)
//{
//  int i,j,k,m,n;
//
//  md->dim=dim;
//  md->cls=cls;
//  md->exist=exist;
//
//  md->mean=(double *)calloc(dim,sizeof(double));
//  md->sigma = matrix_2d_double(dim,dim);
//  md->sigma_inv = matrix_2d_double(dim, dim);
//  return 0;
//}
//
//int cpgauss(GaussModel *md1, GaussModel *md2)
//{
//  int i,j,k,m,n;
//
//  md2->dim=md1->dim;
//  md2->cls=md1->cls;
//  md2->exist=md1->exist;
//  md2->sigma_det=md1->sigma_det;
//
//  for (i=0;i<md1->dim;i++)
//    md2->mean[i]=md1->mean[i];
//  
//  for (i=0;i<md1->dim;i++){
//    for (j=0;j<md1->dim;j++) {
//      md2->sigma[i][j]=md1->sigma[i][j];
//      md2->sigma_inv[i][j]=md1->sigma_inv[i][j];
//    }
//  }
//  return 0;
//}
//
//
//void newhmm(HmmModel *md, int dim, int numst, int numcls, std::vector<int>& stcls)
//{
//  int i,j,k,m,n;
//  int cls, exist=1;
//
//  md->dim=dim;
//  md->numst=numst;
//  md->numcls=numcls;
//  
//  md->stcls=(int *)calloc(numst,sizeof(int));
//  if (numcls==numst || stcls.size()==0) {
//    for (i=0;i<numst;i++) md->stcls[i]=i;
//  }
//  else {
//    for (i=0;i<numst;i++) md->stcls[i]=stcls[i];
//  }
//
//  md->stpdf=(GaussModel **)calloc(numst, sizeof(GaussModel *));
//  for (i=0; i<numst; i++) {
//    md->stpdf[i]=(GaussModel *)calloc(1,sizeof(GaussModel));
//    newgauss(md->stpdf[i], dim, md->stcls[i], exist);
//  }
//  
//  md->a = matrix_2d_double(numst,numst);
//  md->a00=(double *)calloc(numst, sizeof(double));
//}
//
//
//void freehmm(HmmModel **md_pt)
//{
//  int i,j,k,m,n;
//  int cls, exist=1,numst;
//  HmmModel *md;
//
//  md= *md_pt;
//  numst=md->numst;
//
//  for (i=0; i<numst; i++) {
//    free(md->stpdf[i]->mean);
//    free_matrix_2d_double(md->stpdf[i]->sigma, md->dim);
//    free_matrix_2d_double(md->stpdf[i]->sigma_inv, md->dim);
//    free(md->stpdf[i]);
//  }
//  free(md->stpdf);
//
//  free(md->a00);
//  free_matrix_2d_double(md->a, numst);
//
//  free(md->stcls);
//
//  free(md);
//  *md_pt=NULL;
//}
//
//
//void cphmm(HmmModel *md1, HmmModel *md2)
//{
//  int i,j,k,m,n;
//  int cls, exist=1;
//  int numst, numcls,dim;
//
//  md2->dim=dim=md1->dim;
//  md2->numst=numst=md1->numst;
//  md2->numcls=numcls=md1->numcls;
//  
//  for (i=0;i<numst;i++) md2->stcls[i]=md1->stcls[i];
//
//  for (i=0; i<numst; i++) {
//    cpgauss(md1->stpdf[i], md2->stpdf[i]);
//  }
//
//  for (i=0;i<numst;i++) md2->a00[i]=md1->a00[i];
//  for (i=0;i<numst;i++)
//    for (j=0;j<numst;j++)
//      md2->a[i][j]=md1->a[i][j];
//}


double gauss_pdf_log(std::vector<float> &ft, GaussModel &gm, int baseidx)
{
  double res, tpdb, tpdb2;
  int i,j,k,m,n;

  std::vector<double> db_array(gm.dim, 0.0);
  std::vector<double> dif(gm.dim, 0.0);

  m=gm.dim;
  for (int i=0; i<m; i++) {
    dif[i] = ft[baseidx + i] - gm.mean[i];
  }

//  ptrdb1 = db_array;
  for (int i=0; i<m; i++) {
    for (int j=0; j<m; j++) {
      db_array[i] += gm.sigma_inv[i][j] * dif[j];
    }
  }

  tpdb2 = 0.0;
  for (int i=0; i<m; i++){
    tpdb2 += db_array[i]*dif[i];
  }

  tpdb = -((double)(gm.dim))/2.0*LOG_2_PI-0.5*log(gm.sigma_det);
  res = tpdb +(-0.5)*tpdb2;

  return(res);
}

double gauss_pdf(std::vector<float> &ft, GaussModel &gm, int baseidx)
{
  return(exp(gauss_pdf_log(ft, gm, baseidx)));
}

double mix_gauss_pdf_log(std::vector<float> &ft, std::vector<GaussModel> &gmlist, std::vector<double> &prior,
			 int ncmp, int baseidx)
{
  double res, v1,v2;
  int i,j,k,m,n;

  std::vector<double> h(ncmp, 0.0);
  for (i=0;i<ncmp;i++)
    h[i]=gauss_pdf_log(ft, gmlist[i], baseidx);

  v1=h[0];
  for (i=1;i<ncmp;i++) if (h[i]>v1) v1=h[i];

  for (i=0,v2=0.0;i<ncmp;i++) {
    v2 += prior[i]*exp(h[i]-v1);
  }

  if (v2>0.0)  res=v1+log(v2); else res=-HUGE;

  return(res);
}

