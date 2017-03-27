#include "hmm.h"
#include "glog/logging.h"
#include <assert.h>

// see hmm1.in for input format example. specifications will be added later
void HmmModel::read_model(std::string filename){
  FILE *fp = NULL;
  fp = fopen(filename.c_str(), "r");
  assert(fp);
  int dim, numst, tmp;
  char comments[20];
  fscanf(fp, "%d", &dim);
  fscanf(fp, "%d", &numst);
//newhmm(md, dim, numst, numst, NULL);
  this->dim=dim;
  this->numst=numst;
  //we assume numcls==numst, which means each state has 1 Gaussian.
  this->numcls=numst;
  stcls.resize(numst);
  for (int i=0; i<numst; i++)
    stcls[i]=i;
//  stpdf.resize(numst);
  for (int i=0; i<numst; i++) {
    GaussModel tmp(dim, stcls[i], 1);
    stpdf.push_back(tmp);
//    stpdf[i] = GaussModel(dim, stcls[i], 1);
  }
  a.resize(numst);
  for (int i=0; i<numst; i++) {
    a[i].resize(numst);
  }
  a00.resize(numst);
  //skip comments
  fscanf(fp, "%20s", comments);
  for (int i=0; i<numst; i++) {
    fscanf(fp, "%lf ", &this->a00[i]);
  }
  //skip comments
  fscanf(fp, "%20s", comments);
  //read transition matrix a:
  for (int i=0; i<numst; i++) {
    for (int j=0; j<numst; j++)
      fscanf(fp, "%lf ", &this->a[i][j]);
  }
  //read means:
  //skip comments
  fscanf(fp, "%20s", comments);
  for (int i=0; i<numst; i++) {
    for (int j=0; j<dim; j++)
      fscanf(fp, "%lf ", &this->stpdf[i].mean[j]);
  }
  //read sigmas:
  for (int i=0; i<numst; i++) {
    //skip comments
    fscanf(fp, "%20s", comments);
    for (int m=0; m<this->stpdf[i].dim; m++) {
      for (int n=0; n<this->stpdf[i].dim; n++)
        fscanf(fp, "%lf ", &this->stpdf[i].sigma[m][n]);
    }
    tmp=mat_det_inv(this->stpdf[i].sigma, this->stpdf[i].sigma_inv,
                          this->stpdf[i].sigma_det,this->stpdf[i].dim);
//    LOG(INFO)<<"mat_det_inv return value = "<<tmp;
  }
  fclose(fp);
}

void HmmModel::write_model(std::string filename) const{
  FILE *fp = NULL;
  fp = filename.size()>0? fopen(filename.c_str(), "w+") : stdout;
  assert(fp);
  int dim = this->dim, numst = this->numst;
  fprintf(fp, "%d\n", dim);
  fprintf(fp, "%d\n", numst);
  fprintf(fp, "//a00\n");
  for (int i=0; i<numst; i++) {
    fprintf(fp, "%lf ", this->a00[i]);
  }
  fprintf(fp, "\n");
  //write transition matrix a:
  fprintf(fp, "//a\n");
  for (int i=0; i<numst; i++) {
    for (int j=0; j<numst; j++)
      fprintf(fp, "%lf ", this->a[i][j]);
    fprintf(fp, "\n");
  }
  //write means:
  fprintf(fp, "//mus\n");
  for (int i=0; i<numst; i++) {
    for (int j=0; j<dim; j++)
      fprintf(fp, "%lf ", this->stpdf[i].mean[j]);
    fprintf(fp, "\n");
  }
  for (int i=0; i<numst; i++) {
    //skip comments
    fprintf(fp, "//sigmas\n");
    for (int m=0; m<this->stpdf[i].dim; m++) {
      for (int n=0; n<this->stpdf[i].dim; n++)
        fprintf(fp, "%lf ", this->stpdf[i].sigma[m][n]);
      fprintf(fp, "\n");
    }
  }
}



/*------------------------------------------------------------*/
void HmmModel::print_model(std::string filename) const
{
  //Using "" as input filename if you want to print to stdout.
  FILE *outfile = filename.size()>0? fopen(filename.c_str(), "w+") : stdout;
  assert(outfile);
  int dim,numst,numcls;

  dim=this->dim;
  numst=this->numst;
  numcls=this->numcls;
  
  fprintf(outfile, "dim=%d\n", dim);
  fprintf(outfile, "numst=%d\n", numst);
  fprintf(outfile, "numcls=%d\n", numcls);

  fprintf(outfile, "\nState class membership:\n");
  for (int i=0; i<numst; i++)
    fprintf(outfile, "%d ", this->stcls[i]);
  fprintf(outfile, "\n");

  fprintf(outfile, "\nTransition probability a00:\n");
  for (int i=0; i<numst; i++)
    fprintf(outfile, "%8e ", this->a00[i]);
  fprintf(outfile, "\n");

  fprintf(outfile, "Transition probability a:\n");
  for (int i=0; i<numst; i++) {
    for (int j=0; j<numst; j++)
      fprintf(outfile, "%8e ", this->a[i][j]);
    fprintf(outfile, "\n");
  }

  fprintf(outfile, "\nThe Gaussian distributions of states:\n");
  for (int i=0; i<numst; i++) {
    fprintf(outfile, "\nState %d =============================\n", i);
    fprintf(outfile, "exist=%d, dim=%d\n", this->stpdf[i].exist,
	    this->stpdf[i].dim);

    fprintf(outfile, "Mean vector:\n");
    for (int j=0; j<dim; j++)
      fprintf(outfile, "%.5e ", this->stpdf[i].mean[j]);
    fprintf(outfile, "\n");

    fprintf(outfile, "Sigma_det=%e\n", this->stpdf[i].sigma_det);

    fprintf(outfile, "Covariance matrix Sigma:\n");
 
    for (int m=0; m<this->stpdf[i].dim; m++) {
      for (int n=0; n<this->stpdf[i].dim; n++)
        fprintf(outfile, "%.5e ", this->stpdf[i].sigma[m][n]);
      fprintf(outfile, "\n");
    }

    fprintf(outfile, "Covariance matrix inverse Sigma_inv:\n");
 
    for (int m=0; m<this->stpdf[i].dim; m++) {
      for (int n=0; n<this->stpdf[i].dim; n++)
        fprintf(outfile, "%.5e ", this->stpdf[i].sigma_inv[m][n]);
      fprintf(outfile, "\n");
    }
  }
  
}

