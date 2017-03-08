#include "hmm.h"

void forward(float *u, int ncols, double *thetalog, HmmModel *md,
	     double *loglikehd)
// ncols is the length of the sequence
// u[sequence length*dim] stores the sequence of vectors as one array
// thetalog is the log of the forward probabilities
// md is the input HMM
// *loglikehd is the log likelihood for the entire sequence
{
  int i,j,k,l,m,n,ii,jj,mm,nn;
  int precls, curcls;
  int numst,dim;
  double *dbpt,v1,v2,v3,maxv,*buf, *a, **a2;

  numst=md->numst;
  dim=md->dim;

  buf=(double *)calloc(numst,sizeof(double));

  /* treat the first postion */
  a=md->a00;

  for (l=0; l<numst; l++) {
    if (a[l]>0.0)
      thetalog[l]=log(a[l])+gauss_pdf_log(u, md->stpdf[l]);
    else
      thetalog[l]=-HUGE;
  }

  /* treat the rest columns */
  for (jj=1; jj<ncols; jj++) {
    for (l=0;l<numst;l++) {
      buf[l]=thetalog[(jj-1)*numst+l];
    }
    maxv=buf[0];
    for (l=0; l<numst; l++)
      if (buf[l]>maxv) maxv=buf[l];

    a2=md->a;

    for (m=0, mm=jj*numst; m<numst; m++,mm++) {
	v3=gauss_pdf_log(u+jj*dim, md->stpdf[m]);

	for (l=0,v1=0.0;l<numst;l++) {
	  v1+=exp(buf[l]-maxv)*a2[l][m];
	}
	if (v1>0.0)
	  thetalog[mm]=maxv+log(v1)+v3;
	else
	  thetalog[mm]=-HUGE;
    }
  }


  /** compute loglikelihood for the image **/
  v3=0.0;
  dbpt=thetalog+(ncols-1)*numst;
  v3=dbpt[0];
  for (m=0; m<numst; m++) {
    if (dbpt[m]>v3) v3=dbpt[m];
  }
  for (m=0,v1=0.0; m<numst; m++) {
    v1+=exp(dbpt[m]-v3);
  }
  v3+=log(v1);

  *loglikehd=v3;

  free(buf);
}


void backward(float *u, int ncols, double *betalog, HmmModel *md)
// ncols is the length of the sequence
// u[sequence length*dim] stores the sequence of vectors as one array
// betalog is the log of the forward probabilities
// md is the input HMM
{
  int i,j,k,l,m,n,ii,jj,mm,nn;
  int nextcls, curcls;
  int numst,dim;
  double *dbpt,v1,v2,v3,maxv,*buf, *a, **a2;

  numst=md->numst;
  dim=md->dim;

  buf=(double *)calloc(numst,sizeof(double));

  /* treat the last column */
  for (l=0; l<numst; l++) {
    betalog[(ncols-1)*numst+l]=0.0;
  }
  /* treat the rest columns */

  for (jj=ncols-2; jj>=0; jj--) {
    for (l=0; l<numst; l++) {
      buf[l]=betalog[(jj+1)*numst+l]+
	gauss_pdf_log(u+(jj+1)*dim,md->stpdf[l]);
    }

    maxv=buf[0];
    for (l=0; l<numst; l++)
      if (buf[l]>maxv) maxv=buf[l];

    a2=md->a;

    for (m=0, mm=jj*numst; m<numst; m++,mm++) {
      for (l=0,v1=0.0;l<numst;l++) {
	v1+=exp(buf[l]-maxv)*a2[m][l];
      }
      if (v1>0.0)
	betalog[mm]=maxv+log(v1);
      else
	betalog[mm]=-HUGE;
    }
  }

  free(buf);
}


void CompLm(float *u, int ncols, double *thetalog, double *betalog,
	    double *Lm, HmmModel *md)
     /* Lm=double[ncols*numst], space allocated */
{
  int t,i,j,k,l,m,n;
  double *curLm,v1,v2;
  int curcls;
  int numst;

  numst=md->numst;

  for (j=0; j<ncols; j++) {
    curLm=Lm+j*numst;
    for (m=0; m<numst;m++)
      curLm[m]=thetalog[j*numst+m]+betalog[j*numst+m];

    v1=curLm[0];
    for (m=0; m<numst; m++)
      if (curLm[m]>v1) v1=curLm[m];

    for (m=0,v2=0.0;m<numst;m++){
      curLm[m]=exp(curLm[m]-v1);
      v2+=curLm[m];
    }

    for (m=0;m<numst;m++)
      curLm[m]/=v2;
  }
}


void CompHml(float *u, int ncols, double *thetalog, double *betalog,
	      double *Hml, int west, HmmModel *md)
     /* Hml=double[ncols*numst], space allocated */
{
  int t,i,j,k,l,m,n,ll;
  double *curHml,v1,v2;
  int dim,numst;
  double loglikehd, *dbpt;

  dim=md->dim;
  numst=md->numst;

  /* Hml, m is fixed state of the previous neighbor */

  /* compute the log likelihood for the whole sequence */
  loglikehd=0.0;
  dbpt=thetalog+(ncols-1)*numst;

  loglikehd=dbpt[0];
  for (m=0; m<numst; m++) {
      if (dbpt[m]>loglikehd) loglikehd=dbpt[m];
  }

  for (m=0,v1=0.0; m<numst; m++) {
      v1+=exp(dbpt[m]-loglikehd);
  }
  loglikehd+=log(v1);

  for (j=0; j<ncols; j++) {
    curHml=Hml+j*numst;
    if (j==0) {
      for (l=0; l<numst;l++) curHml[l]=0.0;
      continue;
    }

    for (l=0;l<numst;l++) {
      curHml[l]=-loglikehd+thetalog[(j-1)*numst+west]+betalog[j*numst+l]+
	gauss_pdf_log(u+j*dim,md->stpdf[l]);
      curHml[l]=exp(curHml[l])*md->a[west][l];
    }
  }
}

/*-------------------------------*/
/* Viterbi for a single sequence */
/*-------------------------------*/
void viterbi(HmmModel *md, float *u, int len, int *optst, double *inita,
	     double *lastmerit)
// optst stores the optimal sequence of states with maximum posterior
{
  int i,j,k,l,m,n;
  int numst, dim;
  double *merit, **a, *astart;
  int *prest;
  double v1,v2,v3;

  numst=md->numst;
  dim=md->dim;

  prest=(int *)calloc(len*numst,sizeof(int));
  merit=(double *)calloc(len*numst,sizeof(double));

  if (inita==NULL)
    astart=md->a00;
  else
    astart=inita;

  /* treat the first location */
  for (l=0; l<numst; l++) {
    if (astart[l]>0.0) {
      merit[l]=log(astart[l])+gauss_pdf_log(u, md->stpdf[l]);
    }
    else {
      merit[l]=-HUGE;
    }
  }

  /* treat the rest locations */
  for (j=1; j<len; j++) {
    a=md->a;

    for (l=0; l<numst; l++) {
      v1=gauss_pdf_log(u+j*dim, md->stpdf[l]);
      v2=(a[0][l]>0.0)?(merit[(j-1)*numst]+log(a[0][l])):(-HUGE);
      prest[j*numst+l]=0;
      for (m=1; m<numst; m++) {
	v3=(a[m][l]>0.0)?(merit[(j-1)*numst+m]+log(a[m][l])):(-HUGE);
	if (v2<v3) {
	  v2=v3;
	  prest[j*numst+l]=m;
	}
      }
      merit[j*numst+l]=v2+v1;
    }
  }

  m=0;
  v1=merit[(len-1)*numst];
  for (l=1;l<numst;l++) {
    if (merit[(len-1)*numst+l]>v1) {
      v1=merit[(len-1)*numst+l];
      m=l;
    }
  }

  if (lastmerit!=NULL) {
    for (l=0;l<numst;l++) lastmerit[l]=merit[(len-1)*numst+l];
  }

  optst[len-1]=m;
  for (j=len-2; j>=0; j--) {
    optst[j]=prest[(j+1)*numst+optst[j+1]];
  }

  free(prest);
  free(merit);
}


/*---------------------------------------------------------------------*/
/* Viterbi for a single sequence, each component is a gaussian mixture */
/*---------------------------------------------------------------------*/

void formmix(HmmModel *md, double *inita, double **a, double *astart,
	     GaussModel ***pdflist, double **prior, int *nstpercls)
{
  int i,j,k,m,n;
  int numst,numcls;
  int *cls2st;
  double v1,v2,v3;

  numst=md->numst;
  numcls=md->numcls;
  cls2st=(int *)calloc(numcls,sizeof(int));

  for (i=0;i<numcls;i++) {
    for (j=0;j<numst;j++) {
      if (md->stcls[j]==i) {
	cls2st[i]=j;
	break;
      }
    }
  }

  if (inita==NULL) {
    for (i=0;i<numcls;i++) astart[i]=0.0;
    for (j=0;j<numst;j++) astart[md->stcls[j]]+=md->a00[j];
  }
  else {
    for (i=0;i<numcls;i++) astart[i]=0.0;
    for (j=0;j<numst;j++) astart[md->stcls[j]]+=inita[j];
  }

  for (i=0;i<numcls;i++) {
    for (j=0;j<numcls;j++) {
      a[i][j]=0.0;
      for (k=0;k<numst;k++) {
	if (md->stcls[k]==j) {
	  a[i][j] += md->a[cls2st[i]][k];
	}
      }
    }
  }

  for (i=0;i<numcls;i++) nstpercls[i]=0;
  for (j=0;j<numst;j++) nstpercls[md->stcls[j]]++;

  for (i=0;i<numcls;i++) {
    for (j=0,k=0,v1=0.0;j<numst;j++) {
      if (md->stcls[j]==i) {
	prior[i][k]=md->a00[j];
	v1+=prior[i][k];

	pdflist[i][k]=md->stpdf[j];

	k++;
      }
    }
    if (v1>0.0) {
      for (j=0;j<k;j++)	prior[i][j]/=v1;
    }
    else {
      for (j=0;j<k;j++)	prior[i][j]=1.0/(double)k;
    }
  }

  free(cls2st);
}


void viterbicls(HmmModel *md, float *u, int len, int *optst, double *inita,
		double *lastmerit, int *bestnext)
{
  int i,j,k,l,m,n;
  int numst, dim, numcls;
  double *merit, **a, *astart, **prior;
  int *prest;
  double v1,v2,v3;
  GaussModel ***pdflist;
  int *nstpercls;

  numst=md->numst;
  numcls=md->numcls;
  dim=md->dim;

  prest=(int *)calloc(len*numcls,sizeof(int));
  merit=(double *)calloc(len*numcls,sizeof(double));
  pdflist=(GaussModel ***)calloc(numcls,sizeof(GaussModel **));
  for (i=0;i<numcls;i++)
    pdflist[i]=(GaussModel **)calloc(numst,sizeof(GaussModel *));
  a = matrix_2d_double(numcls,numcls);
  astart=(double *)calloc(numcls,sizeof(double));
  prior = matrix_2d_double(numcls,numst);
  nstpercls=(int *)calloc(numcls,sizeof(int));

  formmix(md, inita, a, astart, pdflist, prior, nstpercls);

  /* treat the first location */
  for (l=0; l<numcls; l++) {
    if (astart[l]>0.0) {
      merit[l]=log(astart[l])+
	mix_gauss_pdf_log(u, pdflist[l], prior[l],nstpercls[l]);
    }
    else {
      merit[l]=-HUGE;
    }
  }

  /* treat the rest locations */
  for (j=1; j<len; j++) {
    for (l=0; l<numcls; l++) {
      v1=mix_gauss_pdf_log(u+j*dim, pdflist[l], prior[l], nstpercls[l]);
      v2=(a[0][l]>0.0)?(merit[(j-1)*numcls]+log(a[0][l])):(-HUGE);
      prest[j*numcls+l]=0;
      for (m=1; m<numcls; m++) {
	v3=(a[m][l]>0.0)?(merit[(j-1)*numcls+m]+log(a[m][l])):(-HUGE);
	if (v2<v3) {
	  v2=v3;
	  prest[j*numcls+l]=m;
	}
      }
      merit[j*numcls+l]=v2+v1;
    }
  }

  m=0;
  v1=merit[(len-1)*numcls];
  for (l=1;l<numcls;l++) {
    if (merit[(len-1)*numcls+l]>v1) {
      v1=merit[(len-1)*numcls+l];
      m=l;
    }
  }

  optst[len-1]=m;
  for (j=len-2; j>=0; j--) {
    optst[j]=prest[(j+1)*numcls+optst[j+1]];
  }

  if (lastmerit!=NULL) {
    for (l=0;l<numcls;l++) lastmerit[l]=merit[(len-1)*numcls+l];

    k=-1;
    v1=0.0;
    for (l=0;l<numcls;l++) {
      v2=(a[0][l]>0.0)?(lastmerit[0]+log(a[0][l])):(-HUGE);
      for (m=1; m<numcls; m++) {
	v3=(a[m][l]>0.0)?(lastmerit[m]+log(a[m][l])):(-HUGE);
	if (v2<v3) {
	  v2=v3;
	}
      }
      if (k<0) {
	k=l;
	v1=v2;
      }
      else {
	if (v2>v1) {
	  v1=v2;
	  k=l;
	}
      }
    }

    *bestnext=k;
  }


  free(prest);
  free(merit);
  free_matrix_2d_double(a, numcls);
  free(astart);
  free_matrix_2d_double(prior, numcls);
  free(nstpercls);
  for (i=0;i<numcls;i++) free(pdflist[i]);
  free(pdflist);
}


/*--------------------------------*/
/* Viterbi for multiple sequences */
/*--------------------------------*/

void viterbi_mulseq(HmmModel *md, float **u, int nseq, int *len,
		    int **st)
     /* st=int[nseq][len[?]], space allocated */
{
  int i,j,k,l,m,n;
  int dim;

  dim=md->dim;

  /* treat the rest rows */
  for (i=0; i<nseq; i++) {
    viterbi(md, u[i], len[i], st[i], NULL, NULL);
  }
}


void updatepar_adder(float *u, int ncols, double *thetalog, double *betalog,
		     double loglikehd, HmmModel *md, double *musum,
		     double **mom2sum, double **asum, double *lsum)
{
  int i,j,k,l,m,n,ii,jj,mm,nn,ll,t,k1,k2;
  double *Lm, *Hml;
  double v1,v2,tpdb,epsilon=1.0e-2;
  int dim,numst;

  dim=md->dim;
  numst=md->numst;

  /** allocate space **/
  Hml=(double *)calloc(ncols*numst,sizeof(double));
  Lm=(double *)calloc(ncols*numst,sizeof(double));

  /* Initialization */
  for (i=0; i<numst; i++) {
    lsum[i]=0.0;
    for (j=0; j<dim; j++)
      musum[i*dim+j]=0.0;
    for (j=0; j<dim; j++)
      for (k=0; k<dim; k++)
	mom2sum[i*dim+j][k]=0.0;
  }
  for (j=0; j<numst; j++)
    for (k=0; k<numst; k++)
      asum[j][k]=0.0;

  CompLm(u, ncols, thetalog, betalog, Lm, md);

  for (i=0; i<ncols; i++) {
    for (m=0; m<numst; m++) {
      lsum[m] += Lm[i*numst+m];
      for (jj=0; jj<dim; jj++)
	musum[m*dim+jj]+= Lm[i*numst+m]*u[i*dim+jj];

      /* mom2sum is the second order moment */
      /* covariance is second moment minus the square of the mean */
      for (ii=0; ii<dim; ii++)
	for (jj=0; jj<dim; jj++)
	  mom2sum[m*dim+ii][jj]+=Lm[i*numst+m]*u[i*dim+ii]*u[i*dim+jj];
    }
  }

  for (jj=0; jj<numst; jj++) {
    CompHml(u, ncols, thetalog, betalog, Hml, jj, md);
    for (i=0; i<ncols; i++) {
      for (m=0; m<numst; m++)
	asum[jj][m] += Hml[i*numst+m];
    }
  }

  free(Hml);
  free(Lm);
}


/*-------------------------------*/
/** initialization using kmeans **/
/*-------------------------------*/

void initialize(float **u, int nseq, int *len, int dim, HmmModel *md, int ranflag)
{
  int i,j,k,l,m,n,ii,jj, mm,nn, k1, k2;
  int numst, numdata;
  float *buf,*cdbk;
  int *code;
  double tpdb, epsilon=1.0e-2, lambda=0.5;
  double **sigma, **sigma_inv, *mu, **sigcom;

  /** use kmeans to decide the initial states **/
  numst=md->numst;
  for (i=0,numdata=0;i<nseq;i++) numdata+=len[i];

  buf = (float *)calloc(numdata*dim,sizeof(float));
  code=(int *)calloc(numdata,sizeof(int));
  cdbk=(float *)calloc(numdata*dim,sizeof(float));
  sigcom = matrix_2d_double(dim, dim);

  for (i=0,k=0; i<nseq; i++) {
    for (j=0; j<len[i]; j++) {
      for (jj=0; jj<dim; jj++) buf[k*dim+jj]=u[i][j*dim+jj];
      k++;
    }
  }

  if (!ranflag) {
    lloyd(cdbk,dim,numst,buf,numdata,1.0e-4);
    encode(cdbk,dim,numst,buf,code,numdata);
  }
  else {
    srand48(0);
    for (i=0;i<numdata;i++) {
      code[i]=(int)(drand48()*numst);
      if (code[i]>=numst) code[i]=numst-1;
    }
  }

  /* compute gaussian mean */
  for (m=0; m<numst; m++) {
    mu=md->stpdf[m]->mean;
    for (jj=0; jj<dim; jj++) mu[jj]=0.0;
    for (i=0,n=0; i<numdata; i++) {
      if (code[i]==m) {
        n++;
          for (jj=0; jj<dim; jj++)
            mu[jj]+=buf[i*dim+jj];
      }
    }

    if (n==0) { //don't divide for empty cell
      for (jj=0; jj<dim; jj++)
        mu[jj]=0.0;
    }
    else {
      for (jj=0; jj<dim; jj++)
        mu[jj]/=(double)n;
    }
  }


  // Compute common covariance matrix
  for (k1=0; k1<dim; k1++)
    for (k2=0; k2<dim; k2++)
      sigcom[k1][k2]=0.0;
  for (i=0; i<numdata; i++) {
    for (k1=0; k1<dim; k1++)
      for (k2=k1; k2<dim; k2++)
	sigcom[k1][k2]+=(buf[i*dim+k1]-md->stpdf[code[i]]->mean[k1])*
	  (buf[i*dim+k2]-md->stpdf[code[i]]->mean[k2]);
  }
  for (k1=0; k1<dim; k1++)
    for (k2=k1; k2<dim; k2++) {
      sigcom[k1][k2]/=(double)numdata;
      if (k1!=k2) sigcom[k1][k2]=0.0;  //remove if no spherification

      if (k2!=k1) sigcom[k2][k1]=sigcom[k1][k2];
    }

  /* compute gaussian mean and covariance matrix */
  for (m=0; m<numst; m++) {
    mu=md->stpdf[m]->mean;
    sigma=md->stpdf[m]->sigma;
    sigma_inv=md->stpdf[m]->sigma_inv;

    for (i=0,n=0; i<numdata; i++) {
      if (code[i]==m) 	n++;
    }

    if (n==0) { //don't divide for empty cell
      for (k1=0; k1<dim; k1++)
	for (k2=0; k2<dim; k2++)
	  sigma[k1][k2]=sigcom[k1][k2];  // use common covariance
    }
    else {
      for (k1=0; k1<dim; k1++)
	for (k2=0; k2<dim; k2++)
	  sigma[k1][k2]=0.0;
      for (i=0; i<numdata; i++) {
	if (code[i]==m) {
	  for (k1=0; k1<dim; k1++)
	    for (k2=k1; k2<dim; k2++)
	      sigma[k1][k2]+=(buf[i*dim+k1]-mu[k1])*(buf[i*dim+k2]-mu[k2]);
	}
      }
      for (k1=0; k1<dim; k1++)
	for (k2=k1; k2<dim; k2++) {
	  sigma[k1][k2]/=(double)n;
	  if (k2!=k1) sigma[k2][k1]=sigma[k1][k2];
	}
      for (k1=0; k1<dim; k1++)
	for (k2=0; k2<dim; k2++) {
	  if (!ranflag) //dilate slightly
	    sigma[k1][k2]=sigma[k1][k2]*lambda+(1.1-lambda)*sigcom[k1][k2];
	  else
	    sigma[k1][k2]=sigma[k1][k2]*lambda+(1.0-lambda)*sigcom[k1][k2];

	}
    }

    mm=mat_det_inv_double(sigma,sigma_inv,
			  &(md->stpdf[m]->sigma_det), dim);
    if (mm==2) { /* singular matrix */
      for (k1=0, tpdb=0.0; k1<dim; k1++) tpdb+=sigma[k1][k1];
      tpdb=(tpdb>0.0)?(tpdb/(double)dim*epsilon):epsilon;
      /* modify sigma by adding a scalar matrix */
      for (k1=0; k1<dim; k1++) sigma[k1][k1]+=tpdb;
      mat_det_inv_double(sigma, sigma_inv, &(md->stpdf[m]->sigma_det),dim);
    }
  }

  /** Set the transition probabilities to uniform **/
  tpdb=1.0/(double)numst;
  for (i=0;i<numst;i++) md->a00[i]=tpdb;
  for (i=0;i<numst;i++)
    for (j=0; j<numst; j++)
      md->a[i][j]=tpdb;

  free(buf);
  free(cdbk);
  free(code);
  free_matrix_2d_double(sigcom, dim);

}

/*-----------------------------------------------------------*/
/** compute the log likelihood of sequences under HMM        */
/*-----------------------------------------------------------*/
double comploglike(HmmModel *md, float **u, int nseq, int *len, double *wt)
// weight wt can be turned off by setting it to NULL
{
  int i,j,k,m,n;
  double v1, loglikehd;
  double *thetalog;

  for (i=0,m=0;i<nseq;i++) {
    if (m<len[i]) m=len[i];
  }

  thetalog=(double *)calloc(m*md->numst, sizeof(double));

  for (i=0, loglikehd=0.0; i<nseq; i++) {
    forward(u[i], len[i], thetalog,  md, &v1);
    if (wt==NULL) loglikehd+=v1;
    else loglikehd += wt[i]*v1;
  }

  free(thetalog);
  return(loglikehd);
}


/*-----------------------------------------------------------*/
/** compute the probability of each position existing under  */
/** a particular state given feature vectors for the entire  */
/** sequence and under a given HMM.                          */
/** This subroutine also returns the log likelihood as the   */
/** function comploglike.                                    */
/*-----------------------------------------------------------*/
double classlikehd(HmmModel *md, float **u, int nseq, int *len,
		   double **cprob, double *wt)
     /* cprob[nseq][len[?]] has been allocated with space */
{
  double loglikehd,v1;
  double *thetalog,*betalog;
  double dbtp;
  int i,j,k,m,n,ii,mm;
  int numst;
  int *c, *stcls;

  numst=md->numst;


  for (i=0,mm=0;i<nseq;i++) {
    if (mm<len[i]) mm=len[i];
  }
  thetalog=(double *)calloc(mm, sizeof(double));
  betalog=(double *)calloc(mm, sizeof(double));

  loglikehd=0.0;
  for (ii=0;ii<nseq;ii++) {
    forward(u[ii], len[ii], thetalog, md, &v1);
    backward(u[ii], len[ii], betalog, md);
    CompLm(u[ii], len[ii], thetalog, betalog, cprob[ii], md);

    if (wt==NULL) loglikehd+= v1;
    else loglikehd+= wt[ii]*v1;

    /* normalization */
    for (j=0; j<len[ii]; j++) {
      for (m=0,dbtp=0.0; m<numst; m++)
	dbtp+=cprob[ii][j*numst+m];
      if (dbtp>0.0) {
	for (m=0; m<numst; m++) cprob[ii][j*numst+m]/=dbtp;
      }
      else {
	for (m=0; m<numst; m++) cprob[ii][j*numst+m]=1.0/(double)numst;
      }
    }
  }

  free(thetalog);
  free(betalog);
  return(loglikehd);
}


void transprob(double **asum, int numst, double **amn, int numcls, double **bnl,
	       int *stcls, double *lsum)
{
  int i,j,k,m,n,l;
  double v1,v2,v3;

  for (m=0;m<numcls;m++) {
    for (n=0;n<numcls;n++) {
      amn[m][n]=0.0;
      for (k=0;k<numst;k++) {
	if (stcls[k]!=m) continue;
	for (l=0;l<numst;l++) {
	  if (stcls[l]==n)
	    amn[m][n]+=asum[k][l];
	}
      }
    }
  }

  for (n=0;n<numcls;n++) {
    for (l=0, v1=0.0, k=0; l<numst;l++) {
      if (stcls[l]!=n) bnl[n][l]=0.0;
      else {
	bnl[n][l]=lsum[l];
	v1 += bnl[n][l];
	k++;
      }
    }
    if (k>0) {
      if (v1>0.0)
	for (l=0;l<numst;l++) if (stcls[l]==n) bnl[n][l]/=v1;
      else
	for (l=0;l<numst;l++) if (stcls[l]==n) bnl[n][l]=1.0/(double)k;
    }
  }


  // normalization
  for (k=0; k<numcls; k++) {
    for (l=0,v1=0.0; l<numcls; l++)
      v1+=amn[k][l];
    if (v1>0.0) {
      for (l=0; l<numcls; l++)
	amn[k][l]/=v1;
    }
    else {
      for (l=0; l<numcls; l++)
	amn[k][l]=1.0/(double)numcls;
    }
  }

  // Adjust asum
  for (i=0;i<numst;i++) {
    for (j=0;j<numst;j++) {
      asum[i][j]=amn[stcls[i]][stcls[j]]*bnl[stcls[j]][j];
    }
  }

  for (k=0; k<numst; k++) {
    for (l=0,v1=0.0; l<numst; l++)
      v1+=asum[k][l];
    if (v1>0.0) {
      for (l=0; l<numst; l++)
	asum[k][l]/=v1;
    }
    else {
      for (l=0; l<numst; l++)
	asum[k][l]=1.0/(double)numst;
    }
  }
}


/*---------------------------------------------------------------*/
/** EM estimation assuming that the initial model is set up     **/
/*---------------------------------------------------------------*/
int baumwelch(float **u, int nseq, int *len, HmmModel *md, double *loglikehd,
	      double *lhsumpt, double epsilon, double *wt)
     /* The only outputs are loglikehd, lhsumpt, and updated md */
// Input wt[nseq] gives a weight to each sequence, normally it contains all 1
{
  int i,j,k,l,m,n,mm,k1,m1,t,ite,minite=3, maxrow, maxcol;
  int dim,numst, numcls;
  double ratio=10.0, epsilon2=5.0e-2;
  double oldlhsum, lhsum;
  double *thetalog, *betalog;
  double *musum, *mu, **mom2sum, **mom2, **sigma;
  double **asum, **a, **amn, **bnl, *lsum, *l1img;
  double **sigcom, lambda=LAMBDA;
  int step;
  GaussModel *curg;
  HmmModel *oldmd;
  double v1,tpdb;
  int twomdflag=0;
  int res=0;

  fprintf(stderr, "Inside baumwelch:\n");

  if (nseq==0) return res;

  dim=md->dim;
  numst=md->numst;
  numcls=md->numcls;

  musum = vector_double( numst*dim);
  mu = vector_double( numst*dim);
  mom2sum = matrix_2d_double( numst*dim, dim);
  mom2 = matrix_2d_double( numst*dim, dim);
  sigma = matrix_2d_double( numst*dim, dim);
  sigcom = matrix_2d_double( dim, dim);
  asum = matrix_2d_double( numst, numst);
  a = matrix_2d_double(numst, numst);
  amn = matrix_2d_double( numcls, numcls);
  bnl = matrix_2d_double( numcls, numst);
  lsum = vector_double( numst);
  l1img = vector_double( numst);

  maxcol=len[0];
  for (i=1; i<nseq; i++) {
    if (len[i]>maxcol) maxcol=len[i];
  }
  thetalog=(double *)calloc(maxcol*md->numst, sizeof(double));
  betalog=(double *)calloc(maxcol*md->numst, sizeof(double));


  ite=0;
  twomdflag=0;
  oldlhsum=HUGE;

  while (ite<minite || twomdflag==0 || ratio>epsilon) {
    /* Initialization */
    for (i=0; i<numst; i++) {
      lsum[i]=0.0;
      for (j=0; j<dim; j++)
	musum[i*dim+j]=0.0;
      for (j=0; j<dim; j++)
	for (k=0; k<dim; k++)
	  mom2sum[i*dim+j][k]=0.0;
    }
    for (j=0; j<numst; j++)
      for (l=0; l<numst; l++)
	  asum[j][l]=0.0;

    for (t=0; t<nseq; t++) {
      forward(u[t],len[t],thetalog,md,loglikehd+t);
      backward(u[t],len[t],betalog, md);
      updatepar_adder(u[t],len[t],thetalog, betalog, loglikehd[t],
		      md, mu, mom2, a, l1img);

      for (i=0; i<numst; i++) {
	lsum[i]+=wt[t]*l1img[i];
	for (j=0; j<dim; j++)
	  musum[i*dim+j]+=wt[t]*mu[i*dim+j];
	for (j=0; j<dim; j++)
	  for (k=0; k<dim; k++)
	    mom2sum[i*dim+j][k]+=wt[t]*mom2[i*dim+j][k];
      }
      for (j=0; j<numst; j++)
	for (l=0; l<numst; l++)
	  asum[j][l]+=wt[t]*a[j][l];

    } // for (t=0; ...)

    /* Normalization */
    for (i=0; i<numst; i++) {
      for (j=0; j<dim; j++)
	musum[i*dim+j]/=lsum[i];
      for (j=0; j<dim; j++)
	for (k=0; k<dim; k++)
	  mom2sum[i*dim+j][k]/=lsum[i];
    }

    for (i=0; i<numst; i++) {
      for (j=0; j<dim; j++)
	for (k=0; k<dim; k++)
	  sigma[i*dim+j][k]=mom2sum[i*dim+j][k]-
	    musum[i*dim+j]*musum[i*dim+k];
    }

    // asum adjustment
    transprob(asum,numst,amn,numcls,bnl,md->stcls,lsum);

    for (t=0,lhsum=0.0;t<nseq;t++) lhsum+=wt[t]*loglikehd[t];

    // Judge whether to quit iteration loop
    if (twomdflag>0) {
      ratio=(lhsum-oldlhsum)/fabs(lhsum);
    }
    else {
      ratio=10.0;
    }

    oldlhsum=lhsum;
    ite++;
    fprintf(stderr, "ite=%d, lhsum=%e\n",ite, lhsum);

    if (ratio <= epsilon && ite>=minite) {
      break;
    }

    // exit if empty state appears
    for (i=0,k=0;i<numst;i++) {
      if (lsum[i]==0.0)
	k++;
    }

    if (k) {
      res=1;
      break;
    }

    /*---------------------------*/
    /** update model parameters **/
    /*---------------------------*/
    twomdflag=1;
    for (i=0,v1=0.0; i<numst; i++)
      v1+=lsum[i];
    if (v1>0.0) {
      for (i=0; i<numst; i++)
	md->a00[i]=lsum[i]/v1;
    }
    else {
      for (i=0; i<numst; i++)
//	md->a00[i]=1.0/(double)lsum[i];
          md->a00[i]=1.0/(double)numst;
    }

    for (j=0; j<numst; j++)
      for (k=0; k<numst; k++)
	md->a[j][k]=asum[j][k];

    for (j=0;j<dim;j++)
      for (k=0;k<dim;k++) sigcom[j][k]=0.0;

    for (i=0;i<numst;i++) {
      for (j=0;j<dim;j++)
	for (k=0;k<dim;k++)
	  sigcom[j][k]+=md->a00[i]*sigma[i*dim+j][k];
    }

    for (i=0; i<numst; i++) {
      curg=md->stpdf[i];
      for (j=0; j<dim; j++)
	curg->mean[j]=musum[i*dim+j];
      for (j=0; j<dim; j++)
	for (k=0; k<dim; k++)
	  curg->sigma[j][k]=sigma[i*dim+j][k]*lambda+(1.0-lambda)*sigcom[j][k];

      /* compute the inverse sigma and the determinant of sigma */
      mm=mat_det_inv_double(curg->sigma, curg->sigma_inv,
			    &(curg->sigma_det),dim);

      if (mm==2) { /* singular matrix */
	for (k1=0, tpdb=0.0; k1<dim; k1++) tpdb+=curg->sigma[k1][k1];
	tpdb=(tpdb>0.0)?(tpdb/(double)dim*epsilon2):epsilon2;

	/* modify sigma by adding a scalar matrix */
	for (k1=0; k1<dim; k1++) curg->sigma[k1][k1]+=tpdb;
	mat_det_inv_double(curg->sigma, curg->sigma_inv,
			   &(curg->sigma_det),dim);
      }
    }
  } // while (ite<minite ...)

  *lhsumpt=lhsum;

  free(musum);
  free(mu);
  free_matrix_2d_double(mom2sum, numst*dim);
  free_matrix_2d_double(mom2, numst*dim);
  free_matrix_2d_double(sigma, numst*dim);
  free_matrix_2d_double(sigcom,dim);
  free_matrix_2d_double(asum, numst);
  free_matrix_2d_double(a, numst);
  free_matrix_2d_double(amn, numcls);
  free_matrix_2d_double(bnl, numcls);
  free(lsum);
  free(l1img);
  free(thetalog);
  free(betalog);

  return(res);
}


void hmmfit(float **u, int nseq, int *len, int dim, HmmModel *md, int numst,
	    int *stcls, double *loglikehd, double *lhsumpt,
	    double epsilon, double *wt)
{
  int t,i,j,k,l,m,n;
  double *thiswt;
  /*** Allocate space for the model ***/

  // The stcls[numst] stores the class label for each state. It's only meaningful
  // for HMM with GMM at each state. If the default HMM with a single Gaussian
  // for each state, stcls[] can be set to NULL
  if (stcls!=NULL) {
    for (i=0,m=0;i<numst;i++) if (m<stcls[i]) m=stcls[i];
    m++;
  } else {m=numst;}
  newhmm(md, dim, numst,m,stcls);

  /*** initialize parameters using the states given by kmeans **/
  initialize(u, nseq, len, dim, md,0);

  /** start em iterative estimation **/
  if (wt==NULL) {
    thiswt=(double *)calloc(nseq, sizeof(double));
    for (i=0;i<nseq;i++) thiswt[i]=1.0;
    baumwelch(u, nseq, len, md, loglikehd, lhsumpt, epsilon, thiswt);
    if (thiswt!=NULL) free(thiswt);
  } else {
    baumwelch(u, nseq, len, md, loglikehd, lhsumpt, epsilon, wt);
  }
}

