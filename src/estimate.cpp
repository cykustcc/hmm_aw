#include "glog/logging.h"

#include "hmm.h"
#include <iostream>

void HmmModel::forward(std::vector<float> &u,
                       int ncols,
                       std::vector<double> &thetalog,
                       double &loglikehd)
// ncols is the length of the sequence
// u[sequence length*dim] stores the sequence of vectors as one array
// thetalog is the log of the forward probabilities
// md is the input HMM
// *loglikehd is the log likelihood for the entire sequence
{
  int precls, curcls;
  double v1,v2,v3,maxv;

  std::vector<double> buf(numst,0.0);

  /* treat the first postion */
  for (int l=0; l<numst; l++) {
    if (a00[l]>0.0)
      thetalog[l]=log(a00[l])+stpdf[l].gauss_pdf_log(u, 0);
    else
      thetalog[l]=-HUGE;
  }

  /* treat the rest columns */
  for (int jj=1; jj<ncols; jj++) {
    for (int l=0;l<numst;l++) {
      buf[l]=thetalog[(jj-1)*numst+l];
    }
    maxv=buf[0];
    for (int l=0; l<numst; l++)
      if (buf[l]>maxv) maxv=buf[l];



    for (int m=0, mm=jj*numst; m<numst; m++,mm++) {
      v3=stpdf[m].gauss_pdf_log(u, jj*dim);
      v1=0.0;
      for (int l=0;l<numst;l++) {
        v1+=exp(buf[l]-maxv)*a[l][m];
      }
      if (v1>0.0)
        thetalog[mm]=maxv+log(v1)+v3;
      else
        thetalog[mm]=-HUGE;
    }
  }


  /** compute loglikelihood for the image **/
  v3=0.0;
//  dbpt=thetalog+(ncols-1)*numst;
//  v3=dbpt[0];
  int baseidx = (ncols-1)*numst;
  v3=thetalog[baseidx];
  for (int m=0; m<numst; m++) {
    if (thetalog[baseidx+m]>v3) v3=thetalog[baseidx+m];
  }
  v1=0.0;
  for (int m=0; m<numst; m++) {
    v1+=exp(thetalog[baseidx+m]-v3);
  }
  v3+=log(v1);

  loglikehd=v3;
}


void HmmModel::backward(std::vector<float> &u,
                        int ncols,
                        std::vector<double> &betalog)
// ncols is the length of the sequence
// u[sequence length*dim] stores the sequence of vectors as one array
// betalog is the log of the forward probabilities
// md is the input HMM
{
  int nextcls, curcls;
  double v1 = 0.0,v3  = 0.0,maxv;

  std::vector<double> buf(numst,0.0);

  /* treat the last column */
  for (int l=0; l<numst; l++) {
    betalog[(ncols-1)*numst+l]=0.0;
  }
  /* treat the rest columns */

  for (int jj=ncols-2; jj>=0; jj--) {
    for (int l=0; l<numst; l++) {
      buf[l]=betalog[(jj+1)*numst+l]+
        stpdf[l].gauss_pdf_log(u,(jj+1)*dim);
    }

    maxv=buf[0];
    for (int l=0; l<numst; l++)
      if (buf[l]>maxv) maxv=buf[l];


    for (int m=0, mm=jj*numst; m<numst; m++,mm++) {
      v1 = 0.0;
      for (int l=0;l<numst;l++) {
        v1+=exp(buf[l]-maxv)*a[m][l];
      }
      if (v1>0.0)
        betalog[mm]=maxv+log(v1);
      else
        betalog[mm]=-HUGE;
    }
  }

}


void HmmModel::CompLm(std::vector<float> &u,
                      int ncols,
                      std::vector<double> &thetalog,
                      std::vector<double> &betalog,
                      std::vector<double> &Lm)
     /* Lm=double[ncols*numst], space allocated */
{
  double v1,v2;

  for (int j=0; j<ncols; j++) {
    int baseidx = j*numst;
    for (int m=0; m<numst;m++)
      Lm[baseidx + m]=thetalog[j*numst+m]+betalog[j*numst+m];

    v1=Lm[baseidx];
    for (int m=0; m<numst; m++)
      if (Lm[baseidx + m]>v1) v1=Lm[baseidx + m];
    
    v2=0.0;
    for (int m=0;m<numst;m++){
      Lm[baseidx + m]=exp(Lm[baseidx + m]-v1);
      v2+=Lm[baseidx + m];
    }

    for (int m=0;m<numst;m++)
      Lm[baseidx+m]/=v2;
  }
}


void HmmModel::CompHml(std::vector<float> &u,
                       int ncols,
                       std::vector<double> &thetalog,
                       std::vector<double> &betalog,
                       std::vector<double> &Hml,
                       int west)
     /* Hml=double[ncols*numst], space allocated */
{
  double v1,v2;
  double loglikehd;
  
  /* Hml, m is fixed state of the previous neighbor */

  /* compute the log likelihood for the whole sequence */
  loglikehd=0.0;
  int baseidx = (ncols-1)*numst;

  loglikehd=thetalog[baseidx];
  for (int m=0; m<numst; m++) {
      if (thetalog[baseidx + m]>loglikehd) loglikehd=thetalog[baseidx + m];
  }
  
  v1=0.0;
  for (int m=0; m<numst; m++) {
      v1+=exp(thetalog[baseidx + m]-loglikehd);
  }
  loglikehd+=log(v1);

  for (int j=0; j<ncols; j++) {
    int baseidx = j*numst;
    if (j==0) {
      for (int l=0; l<numst;l++) Hml[baseidx + l]=0.0;
      continue;
    }

    for (int l=0;l<numst;l++) {
      Hml[baseidx + l]=-loglikehd+thetalog[(j-1)*numst+west]+betalog[j*numst+l]+
        stpdf[l].gauss_pdf_log(u, j*dim);
      Hml[baseidx + l]=exp(Hml[baseidx + l])*a[west][l];
    }
  }
}

/*-------------------------------*/
/* Viterbi for a single sequence */
/*-------------------------------*/
void HmmModel::viterbi(std::vector<float> &u,
                       int len,
                       std::vector<int> &optst,
                       std::vector<double> &inita,
                       std::vector<double> &lastmerit)
// optst stores the optimal sequence of states with maximum posterior
{
  int m;
  double v1,v2,v3;


  std::vector<int> prest(len*numst,0);
  std::vector<double> merit(len*numst,0.0);

  /* treat the first location */
  if (inita.size()==0){
    for (int l=0; l<numst; l++) {
      if (a00[l]>0.0) {
        merit[l]=log(a00[l])+stpdf[l].gauss_pdf_log(u, 0);
      }
      else {
        merit[l]=-HUGE;
      }
    }
  }else{
    for (int l=0; l<numst; l++) {
      if (inita[l]>0.0) {
        merit[l]=log(inita[l])+stpdf[l].gauss_pdf_log(u, 0);
      }
      else {
        merit[l]=-HUGE;
      }
    }
  }

  /* treat the rest locations */
  for (int j=1; j<len; j++) {

    for (int l=0; l<numst; l++) {
      v1=stpdf[l].gauss_pdf_log(u, j*dim);
      v2=(a[0][l]>0.0)?(merit[(j-1)*numst]+log(a[0][l])):(-HUGE);
      prest[j*numst+l]=0;
      for (int m=1; m<numst; m++) {
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
  for (int l=1;l<numst;l++) {
    if (merit[(len-1)*numst+l]>v1) {
      v1=merit[(len-1)*numst+l];
      m=l;
    }
  }

  if (lastmerit.size()!=0) {
    for (int l=0;l<numst;l++) lastmerit[l]=merit[(len-1)*numst+l];
  }

  optst[len-1]=m;
  for (int j=len-2; j>=0; j--) {
    optst[j]=prest[(j+1)*numst+optst[j+1]];
  }
}


/*---------------------------------------------------------------------*/
/* Viterbi for a single sequence, each component is a gaussian mixture */
/*---------------------------------------------------------------------*/

void HmmModel::formmix(std::vector<double> &inita,
                       std::vector<std::vector<double>> &tm,
                       std::vector<double> &astart,
                       std::vector<std::vector<GaussModel>> &pdflist,
                       std::vector<std::vector<double>> &prior,
                       std::vector<int> &nstpercls)
{
  int k;
  double v1;
  
  std::vector<int> cls2st(numcls,0);

  for (int i=0;i<numcls;i++) {
    for (int j=0;j<numst;j++) {
      if (stcls[j]==i) {
        cls2st[i]=j;
        break;
      }
    }
  }

  if (inita.size() == 0) {
    for (int i=0;i<numcls;i++) astart[i]=0.0;
    for (int j=0;j<numst;j++) astart[stcls[j]]+=a00[j];
  }
  else {
    for (int i=0;i<numcls;i++) astart[i]=0.0;
    for (int j=0;j<numst;j++) astart[stcls[j]]+=inita[j];
  }

  for (int i=0;i<numcls;i++) {
    for (int j=0;j<numcls;j++) {
      tm[i][j]=0.0;
      for (int k=0;k<numst;k++) {
        if (stcls[k]==j) {
          tm[i][j] += a[cls2st[i]][k];
        }
      }
    }
  }

  for (int i=0;i<numcls;i++) nstpercls[i]=0;
  for (int j=0;j<numst;j++) nstpercls[stcls[j]]++;

  for (int i=0;i<numcls;i++) {
    v1=0.0;
    for (int j=0,k=0;j<numst;j++) {
      if (stcls[j]==i) {
        prior[i][k]=a00[j];
        v1+=prior[i][k];

        pdflist[i][k]=stpdf[j];

        k++;
      }
    }
    if (v1>0.0) {
      for (int j=0;j<k;j++)	prior[i][j]/=v1;
    }
    else {
      for (int j=0;j<k;j++)	prior[i][j]=1.0/(double)k;
    }
  }
}


void HmmModel::viterbicls(std::vector<float> &u,
                          int len,
                          std::vector<int> &optst,
                          std::vector<double> &inita,
                          std::vector<double> &lastmerit,
                          int &bestnext)
{
  int i,j,k,l,m,n;
  double v1,v2,v3;

  std::vector<int> prest(len*numcls, 0);
  std::vector<double> merit(len*numcls,0.0);
  std::vector<std::vector<GaussModel>> pdflist(numcls, std::vector<GaussModel>(numcls));
  std::vector<std::vector<double>> tm(numcls, std::vector<double>(numcls, 0.0));
  std::vector<double> astart(numcls,0.0);
  std::vector<std::vector<double>> prior(numcls, std::vector<double>(numcls, 0.0));
  std::vector<int> nstpercls(numcls,0);

  formmix(inita, tm, astart, pdflist, prior, nstpercls);

  /* treat the first location */
  for (l=0; l<numcls; l++) {
    if (astart[l]>0.0) {
      merit[l]=log(astart[l])+
      mix_gauss_pdf_log(u, pdflist[l], prior[l],nstpercls[l], 0);
    }
    else {
      merit[l]=-HUGE;
    }
  }

  /* treat the rest locations */
  for (j=1; j<len; j++) {
    for (l=0; l<numcls; l++) {
      v1=mix_gauss_pdf_log(u, pdflist[l], prior[l], nstpercls[l], j*dim);
      v2=(a[0][l]>0.0)?(merit[(j-1)*numcls]+log(a[0][l])):(-HUGE);
      prest[j*numcls+l]=0;
      for (m=1; m<numcls; m++) {
        v3=(a[m][l]>0.0)?(merit[(j-1)*numcls+m]+log(tm[m][l])):(-HUGE);
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

  if (lastmerit.size()!=0) {
    for (l=0;l<numcls;l++) lastmerit[l]=merit[(len-1)*numcls+l];

    k=-1;
    v1=0.0;
    for (l=0;l<numcls;l++) {
      v2=(tm[0][l]>0.0)?(lastmerit[0]+log(tm[0][l])):(-HUGE);
      for (m=1; m<numcls; m++) {
        v3=(tm[m][l]>0.0)?(lastmerit[m]+log(tm[m][l])):(-HUGE);
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

    bestnext=k;
  }

}


/*--------------------------------*/
/* Viterbi for multiple sequences */
/*--------------------------------*/

void HmmModel::viterbi_mulseq(std::vector<std::vector<float>> &u,
                              int nseq,
                              std::vector<int> &len,
                              std::vector<std::vector<int>> &st)
     /* st=int[nseq][len[?]], space allocated */
{
  /* treat the rest rows */
  for (int i=0; i<nseq; i++) {
    std::vector<double> inita;
    std::vector<double> lastmerit;
    viterbi(u[i], len[i], st[i], inita, lastmerit);
  }
}


void HmmModel::updatepar_adder(std::vector<float> &u,
                               int ncols,
                               std::vector<double> &thetalog,
                               std::vector<double> &betalog,
                               double loglikehd,
                               std::vector<double> &musum,
                               std::vector<std::vector<double>> &mom2sum,
                               std::vector<std::vector<double>> &asum,
                               std::vector<double> &lsum)
{
  /** allocate space **/
  std::vector<double> Hml(ncols*numst, 0.0);
  std::vector<double> Lm(ncols*numst, 0.0);

  /* Initialization */
  for (int i=0; i<numst; i++) {
    lsum[i]=0.0;
    for (int j=0; j<dim; j++)
      musum[i*dim+j]=0.0;
    for (int j=0; j<dim; j++)
      for (int k=0; k<dim; k++)
        mom2sum[i*dim+j][k]=0.0;
  }
  for (int j=0; j<numst; j++)
    for (int k=0; k<numst; k++)
      asum[j][k]=0.0;
  
  CompLm(u, ncols, thetalog, betalog, Lm);

  for (int i=0; i<ncols; i++) {
    for (int m=0; m<numst; m++) {
      lsum[m] += Lm[i*numst+m];
      for (int jj=0; jj<dim; jj++)
        musum[m*dim+jj] += Lm[i*numst+m]*u[i*dim+jj];

      /* mom2sum is the second order moment */
      /* covariance is second moment minus the square of the mean */
      for (int ii=0; ii<dim; ii++)
        for (int jj=0; jj<dim; jj++)
          mom2sum[m*dim+ii][jj]+=Lm[i*numst+m]*u[i*dim+ii]*u[i*dim+jj];
    }
  }

  for (int jj=0; jj<numst; jj++) {
    CompHml(u, ncols, thetalog, betalog, Hml, jj);
    for (int i=0; i<ncols; i++) {
      for (int m=0; m<numst; m++)
        asum[jj][m] += Hml[i*numst+m];
    }
  }
}


/*-------------------------------*/
/** initialization using kmeans **/
/*-------------------------------*/

void HmmModel::initialize(std::vector<std::vector<float>> &u,
                          int nseq,
                          std::vector<int> &len,
                          int dim,
                          int ranflag)
{
  int i,j,k,l,m,n,ii,jj, mm,nn, k1, k2;
  int numdata;
  double tpdb, epsilon=1.0e-2, lambda=0.5;

  /** use kmeans to decide the initial states **/
  for (i=0,numdata=0;i<nseq;i++) numdata+=len[i];

  std::vector<float> buf(numdata*dim, 0.0);
  std::vector<int> code(numdata, 0);
  std::vector<float> cdbk(numdata*dim, 0.0);
  std::vector<std::vector<double>> sigcom(dim, std::vector<double>(dim, 0.0));

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
    for (jj=0; jj<dim; jj++)
      stpdf[m].mean[jj]=0.0;
    for (i=0,n=0; i<numdata; i++) {
      if (code[i]==m) {
        n++;
        for (jj=0; jj<dim; jj++)
          stpdf[m].mean[jj]+=buf[i*dim+jj];
      }
    }

    if (n==0) { //don't divide for empty cell
      for (jj=0; jj<dim; jj++)
        stpdf[m].mean[jj] = 0.0;
    }
    else {
      for (jj=0; jj<dim; jj++)
        stpdf[m].mean[jj] /= (double)n;
    }
  }


  // Compute common covariance matrix
  for (k1=0; k1<dim; k1++)
    for (k2=0; k2<dim; k2++)
      sigcom[k1][k2]=0.0;
  for (i=0; i<numdata; i++) {
    for (k1=0; k1<dim; k1++)
      for (k2=k1; k2<dim; k2++)
        sigcom[k1][k2]+=(buf[i*dim+k1]-stpdf[code[i]].mean[k1])*
        (buf[i*dim+k2]-stpdf[code[i]].mean[k2]);
  }

  for (k1=0; k1<dim; k1++)
    for (k2=k1; k2<dim; k2++) {
      sigcom[k1][k2]/=(double)numdata;
      if (k1!=k2) sigcom[k1][k2]=0.0;  //remove if no spherification

      if (k2!=k1) sigcom[k2][k1]=sigcom[k1][k2];
    }

  /* compute gaussian mean and covariance matrix */
  for (m=0; m<numst; m++) {
    for (i=0,n=0; i<numdata; i++) {
      if (code[i]==m) 	n++;
    }

    if (n==0) { //don't divide for empty cell
      for (k1=0; k1<dim; k1++)
        for (k2=0; k2<dim; k2++)
          stpdf[m].sigma[k1][k2]=sigcom[k1][k2];  // use common covariance
    }
    else {
      for (k1=0; k1<dim; k1++)
        for (k2=0; k2<dim; k2++)
          stpdf[m].sigma[k1][k2]=0.0;
      for (i=0; i<numdata; i++) {
        if (code[i]==m) {
          for (k1=0; k1<dim; k1++)
            for (k2=k1; k2<dim; k2++)
              stpdf[m].sigma[k1][k2]+=
              (buf[i*dim+k1]-stpdf[m].mean[k1])*
              (buf[i*dim+k2]-stpdf[m].mean[k2]);
        }
      }
      for (k1=0; k1<dim; k1++)
        for (k2=k1; k2<dim; k2++) {
          stpdf[m].sigma[k1][k2]/=(double)n;
          if (k2!=k1)
            stpdf[m].sigma[k2][k1]=stpdf[m].sigma[k1][k2];
        }
      for (k1=0; k1<dim; k1++)
        for (k2=0; k2<dim; k2++) {
          if (!ranflag) //dilate slightly
            stpdf[m].sigma[k1][k2]=stpdf[m].sigma[k1][k2]*lambda+(1.1-lambda)*sigcom[k1][k2];
          else
            stpdf[m].sigma[k1][k2]=stpdf[m].sigma[k1][k2]*lambda+(1.0-lambda)*sigcom[k1][k2];

        }
    }
    
    mm = mat_det_inv(stpdf[m].sigma,
                     stpdf[m].sigma_inv,
                     stpdf[m].sigma_det,
                     dim);
    if (mm==2) { /* singular matrix */
      for (k1=0, tpdb=0.0; k1<dim; k1++)
        tpdb+=stpdf[m].sigma[k1][k1];
      tpdb=(tpdb>0.0)?(tpdb/(double)dim*epsilon):epsilon;
      /* modify sigma by adding a scalar matrix */
      for (k1=0; k1<dim; k1++)
        stpdf[m].sigma[k1][k1]+=tpdb;
      mat_det_inv(stpdf[m].sigma,
                  stpdf[m].sigma_inv,
                  stpdf[m].sigma_det,
                  dim);
    }
  }

  /** Set the transition probabilities to uniform **/
  tpdb=1.0/(double)numst;
  for (i=0;i<numst;i++)
    a00[i]=tpdb;
  for (i=0;i<numst;i++)
    for (j=0; j<numst; j++)
      a[i][j]=tpdb;
}

/*-----------------------------------------------------------*/
/** compute the log likelihood of sequences under HMM        */
/*-----------------------------------------------------------*/
double HmmModel::comploglike(std::vector<std::vector<float>> &u,
                             int nseq,
                             std::vector<int> &len,
                             std::vector<double>& wt)
// weight wt can be turned off by setting it to NULL
{
  // max length:
  int m = *std::max_element(len.begin(), len.end());
  double v1, loglikehd;
  std::vector<double> thetalog(m*numst, 0.0);

  for (int i=0, loglikehd=0.0; i<nseq; i++) {
    forward(u[i], len[i], thetalog, v1);
    if (wt.size() == 0) loglikehd+=v1;
    else loglikehd += wt[i]*v1;
  }

  return(loglikehd);
}


/*-----------------------------------------------------------*/
/** compute the probability of each position existing under  */
/** a particular state given feature vectors for the entire  */
/** sequence and under a given HMM.                          */
/** This subroutine also returns the log likelihood as the   */
/** function comploglike.                                    */
/*-----------------------------------------------------------*/
double HmmModel::classlikehd(std::vector<std::vector<float>> &u,
                             int nseq,
                             std::vector<int> &len,
                             std::vector<std::vector<double>> &cprob,
                             std::vector<double>& wt)
     /* cprob[nseq][len[?]] has been allocated with space */
{
  double loglikehd,v1;
  double dbtp;
  int i,j,k,m,n,ii,mm;
  int *c, *stcls;



  for (i=0,mm=0;i<nseq;i++) {
    if (mm<len[i]) mm=len[i];
  }
  std::vector<double> thetalog(mm, 0.0);
  std::vector<double> betalog(mm, 0.0);

  loglikehd=0.0;
  for (ii=0;ii<nseq;ii++) {
    forward(u[ii], len[ii], thetalog, v1);
    backward(u[ii], len[ii], betalog);
    CompLm(u[ii], len[ii], thetalog, betalog, cprob[ii]);

    if (wt.size() == 0) loglikehd+= v1;
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
  return(loglikehd);
}


void transprob(std::vector<std::vector<double>> &asum,
               int numst,
               std::vector<std::vector<double>> &amn,
               int numcls,
               std::vector<std::vector<double>> &bnl,
               std::vector<int> &stcls,
               std::vector<double> &lsum)
{
  int k, l;
  double v1;

  for (int m=0;m<numcls;m++) {
    for (int n=0;n<numcls;n++) {
      amn[m][n]=0.0;
      for (int k=0;k<numst;k++) {
        if (stcls[k]!=m) continue;
        for (int l=0;l<numst;l++) {
          if (stcls[l]==n)
            amn[m][n]+=asum[k][l];
        }
      }
    }
  }

  for (int n=0;n<numcls;n++) {
    v1=0.0;
    for (int l=0, k=0; l<numst;l++) {
      if (stcls[l]!=n)
        bnl[n][l]=0.0;
      else {
        bnl[n][l]=lsum[l];
        v1 += bnl[n][l];
        k++;
      }
    }
    if (k>0) {
      if (v1>0.0){
        for (int l=0;l<numst;l++)
          if (stcls[l]==n)
            bnl[n][l]/=v1;
      }else{
          for (int l=0;l<numst;l++)
            if (stcls[l]==n)
              bnl[n][l]=1.0/(double)k;
      }
    }
  }


  // normalization
  for (int k=0; k<numcls; k++) {
    v1=0.0;
    for (int l=0; l<numcls; l++)
      v1+=amn[k][l];
    if (v1>0.0) {
      for (int l=0; l<numcls; l++)
        amn[k][l]/=v1;
    }
    else {
      for (int l=0; l<numcls; l++)
        amn[k][l]=1.0/(double)numcls;
    }
  }

  // Adjust asum
  for (int i=0;i<numst;i++) {
    for (int j=0;j<numst;j++) {
      asum[i][j]=amn[stcls[i]][stcls[j]]*bnl[stcls[j]][j];
    }
  }

  for (int k=0; k<numst; k++) {
    v1=0.0;
    for (int l=0; l<numst; l++)
      v1+=asum[k][l];
    if (v1>0.0) {
      for (int l=0; l<numst; l++)
        asum[k][l]/=v1;
    }
    else {
      for (int l=0; l<numst; l++)
        asum[k][l]=1.0/(double)numst;
    }
  }
}


/*---------------------------------------------------------------*/
/** EM estimation assuming that the initial model is set up     **/
/*---------------------------------------------------------------*/
int HmmModel::baumwelch(std::vector<std::vector<float>> &u,
                        int nseq,
                        std::vector<int> &len,
                        std::vector<double> &loglikehd,
                        double &lhsumpt,
                        double epsilon,
                        std::vector<double> &wt,
                        bool forcediag)
     /* The only outputs are loglikehd, lhsumpt, and updated md */
// Input wt[nseq] gives a weight to each sequence, normally it contains all 1
{
  int i,j,k,l,m,n,mm,k1,m1,t,ite,minite=3, maxrow, maxcol;
  double ratio=10.0, epsilon2=5.0e-2;
  double oldlhsum = 0.0, lhsum = 0.0;
  double lambda=LAMBDA;
  GaussModel *curg;
  double v1,tpdb=0;
  int twomdflag=0;
  int res=0;

  fprintf(stderr, "Inside baumwelch:\n");

  if (nseq==0) return res;

  std::vector<double> musum( numst*dim, 0.0);
  std::vector<double> mu( numst*dim, 0.0);
  std::vector<std::vector<double>> mom2sum( numst*dim, std::vector<double>(dim,0.0));
  std::vector<std::vector<double>> mom2( numst*dim, std::vector<double>(dim,0.0));
  std::vector<std::vector<double>> sigma( numst*dim, std::vector<double>(dim,0.0));
  std::vector<std::vector<double>> sigcom( dim, std::vector<double>(dim,0.0));
  std::vector<std::vector<double>> asum( numst, std::vector<double>(numst,0.0));
  std::vector<std::vector<double>> tm( numst, std::vector<double>(numst,0.0));
  std::vector<std::vector<double>> amn( numst, std::vector<double>(numst,0.0));
  std::vector<std::vector<double>> bnl( numst, std::vector<double>(numst,0.0));
  std::vector<double>lsum( numst, 0.0);
  std::vector<double>l1img( numst, 0.0);

  maxcol=len[0];
  for (int i=1; i<nseq; i++) {
    if (len[i]>maxcol) maxcol=len[i];
  }
  std::vector<double> thetalog(maxcol*numst, 0.0);
  std::vector<double> betalog(maxcol*numst, 0.0);

  ite=0;
  twomdflag=0;
  oldlhsum=HUGE;

  while (ite<minite || twomdflag==0 || ratio>epsilon) {
    /* Initialization */
    for (int i=0; i<numst; i++) {
      lsum[i]=0.0;
      for (int j=0; j<dim; j++)
        musum[i*dim+j]=0.0;
      for (int j=0; j<dim; j++)
        for (int k=0; k<dim; k++)
          mom2sum[i*dim+j][k]=0.0;
    }
    for (int j=0; j<numst; j++)
      for (int l=0; l<numst; l++)
        asum[j][l]=0.0;
    for (int t=0; t<nseq; t++) {
      forward(u[t],len[t],thetalog,loglikehd[t]);
      backward(u[t],len[t],betalog);
      updatepar_adder(u[t],len[t],thetalog, betalog, loglikehd[t],
           mu, mom2, tm, l1img);

      for (int i=0; i<numst; i++) {
        lsum[i]+=wt[t]*l1img[i];
        for (int j=0; j<dim; j++)
          musum[i*dim+j]+=wt[t]*mu[i*dim+j];
        for (int j=0; j<dim; j++)
          for (int k=0; k<dim; k++)
            mom2sum[i*dim+j][k]+=wt[t]*mom2[i*dim+j][k];
      }
      for (int j=0; j<numst; j++)
        for (int l=0; l<numst; l++)
          asum[j][l]+=wt[t]*tm[j][l];

    } // for (t=0; ...)
    
    /* Normalization */
    for (int i=0; i<numst; i++) {
      for (int j=0; j<dim; j++)
        musum[i*dim+j]/=lsum[i];
      for (int j=0; j<dim; j++)
        for (int k=0; k<dim; k++)
        mom2sum[i*dim+j][k]/=lsum[i];
    }
    
    for (int i=0; i<numst; i++) {
      for (int j=0; j<dim; j++){
        for (int k=0; k<dim; k++){
          sigma[i*dim+j][k]=mom2sum[i*dim+j][k]-
          musum[i*dim+j]*musum[i*dim+k];
        }
      }
    }
    // asum adjustment
    transprob(asum,numst,amn,numcls,bnl,stcls,lsum);

    lhsum=0.0;
    for (int t=0;t<nseq;t++) lhsum+=wt[t]*loglikehd[t];

    // Judge whether to quit iteration loop
    if (twomdflag>0) {
      ratio=(lhsum-oldlhsum)/fabs(lhsum);
    }
    else {
      ratio=10.0;
    }

    oldlhsum=lhsum;
    ite++;

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
    v1=0.0;
    for (int i=0; i<numst; i++)
      v1+=lsum[i];
    if (v1>0.0) {
      for (int i=0; i<numst; i++)
        a00[i]=lsum[i]/v1;
    }
    else {
      for (int i=0; i<numst; i++)
      //	md->a00[i]=1.0/(double)lsum[i];
          a00[i]=1.0/(double)numst;
    }

    for (int j=0; j<numst; j++)
      for (int k=0; k<numst; k++)
        a[j][k]=asum[j][k];

    for (int j=0;j<dim;j++)
      for (int k=0;k<dim;k++)
        sigcom[j][k]=0.0;

    for (int i=0;i<numst;i++) {
      for (int j=0;j<dim;j++)
        for (int k=0;k<dim;k++)
          sigcom[j][k] += a00[i]*sigma[i*dim+j][k];
    }

    for (int i=0; i<numst; i++) {
      curg = &stpdf[i];
      for (int j=0; j<dim; j++){
        curg->mean[j]=musum[i*dim+j];
      }
      for (int j=0; j<dim; j++){
        for (int k=0; k<dim; k++){
          curg->sigma[j][k]=sigma[i*dim+j][k]*lambda+(1.0-lambda)*sigcom[j][k];
          if(forcediag && j!=k)
            curg->sigma[j][k] = 0;
        }
      }

      /* compute the inverse sigma and the determinant of sigma */
      mm=mat_det_inv(curg->sigma, curg->sigma_inv,
          curg->sigma_det,dim);

      if (mm==2) { /* singular matrix */
        for (k1=0, tpdb=0.0; k1<dim; k1++)
          tpdb+=curg->sigma[k1][k1];
        tpdb=(tpdb>0.0)?(tpdb/(double)dim*epsilon2):epsilon2;

        /* modify sigma by adding a scalar matrix */
        for (int k1=0; k1<dim; k1++)
          curg->sigma[k1][k1]+=tpdb;
        mat_det_inv(curg->sigma, curg->sigma_inv,
           curg->sigma_det,dim);
      }
    }
  } // while (ite<minite ...)

  lhsumpt=lhsum;
  return(res);
}


void HmmModel::hmmfit(std::vector<std::vector<float>> &u,
                      int nseq,
                      std::vector<int> &len,
                      std::vector<double> &loglikehd,
                      double &lhsumpt,
                      double epsilon,
                      std::vector<double> &wt,
                      bool forcediag)
{
  // The stcls[numst] stores the class label for each state. It's only meaningful
  // for HMM with GMM at each state. If the default HMM with a single Gaussian
  // for each state, stcls[] can be set to NULL
  
  /*** initialize parameters using the states given by kmeans **/
  initialize(u, nseq, len, dim, 0);

  /** start em iterative estimation **/
  if (wt.size() == 0) {
    std::vector<double> thiswt(nseq, 1.0);
    baumwelch(u, nseq, len, loglikehd, lhsumpt, epsilon, thiswt, forcediag);
  } else {
    baumwelch(u, nseq, len, loglikehd, lhsumpt, epsilon, wt, forcediag);
  }
}

