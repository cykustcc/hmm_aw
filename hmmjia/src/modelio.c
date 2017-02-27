#include "hmm.h"

unsigned char write_model(HmmModel *md, FILE *outfile)
{
  int i,j,k,m,n;
  int dim,numst,numcls;

  dim=md->dim;
  numst=md->numst;
  numcls=md->numcls;
  
  if (fwrite(&dim, sizeof(int), 1, outfile)!=1)
    {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }

  if (fwrite(&numst, sizeof(int), 1, outfile)!=1)
    {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }

  if (fwrite(&numcls, sizeof(int), 1, outfile)!=1)
    {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }

  if (fwrite(md->stcls, sizeof(int), numst, outfile)!=numst)
    {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }


  if (fwrite(md->a00, sizeof(double),numst, outfile)!=numst)
    {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }

  for (i=0; i<numst; i++) {
    if (fwrite(md->a[i], sizeof(double),numst, outfile)!=numst)
      {
	fprintf(stderr, "**ERROR writing out data\n");
	return(0);
      }    
  }

  for (i=0; i<numst; i++) {
    if (fwrite(&(md->stpdf[i]->exist), sizeof(int),1,outfile)!=1) {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }

    if (fwrite(&(md->stpdf[i]->dim), sizeof(int), 1, 
	       outfile)!=1) {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    } 
    
    if (fwrite(md->stpdf[i]->mean, sizeof(double), md->stpdf[i]->dim, 
	       outfile)!=md->stpdf[i]->dim) {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }
 
    if (fwrite(&(md->stpdf[i]->sigma_det), sizeof(double),1,
	       outfile)!=1) {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }
 
    for (m=0; m<md->stpdf[i]->dim; m++) {
      if (fwrite(md->stpdf[i]->sigma[m], sizeof(double), md->stpdf[i]->dim,
		 outfile)!=md->stpdf[i]->dim) {
	fprintf(stderr, "**ERROR writing out data\n");
	return(0);
      } 
    }

    for (m=0; m<md->stpdf[i]->dim; m++) {
      if (fwrite(md->stpdf[i]->sigma_inv[m], sizeof(double), 
		 md->stpdf[i]->dim, outfile)!=md->stpdf[i]->dim) {
	fprintf(stderr, "**ERROR writing out data\n");
	return(0);
      } 
    }
  }
  
  return(1);
}

/*------------------------------------------------------------*/

unsigned char read_model(HmmModel *md, FILE *infile)
{
  int i,j,k,m,n;
  int dim,numst,numcls;
  if (fread(&dim, sizeof(int), 1, infile)!=1)
    {
      fprintf(stderr, "**ERROR reading in model dim \n");
      return(0);
    }
  
  if (fread(&numst, sizeof(int), 1, infile)!=1)
    {
      fprintf(stderr, "**ERROR reading in model numst \n");
      return(0);
    }

  if (fread(&numcls, sizeof(int), 1, infile)!=1)
    {
      fprintf(stderr, "**ERROR writing out data numcls \n");
      return(0);
    }

  md->dim=dim;
  md->numst=numst;
  md->numcls=numcls;
  
  md->stcls=(int *)calloc(numst,sizeof(int));
  if (fread(md->stcls, sizeof(int), numst, infile)!=numst)
    {
      fprintf(stderr, "**ERROR writing out data\n");
      return(0);
    }

  md->a00=(double *)calloc(numst,sizeof(double));
  if (md->a00==NULL) {
    fprintf(stderr, "Error allocate space while reading in model\n");
    exit(1);
  }

  if (fread(md->a00, sizeof(double),numst, infile)!=numst)
    {
      fprintf(stderr, "**ERROR reading in model\n");
      return(0);
    }

  md->a=(double **)calloc(numst,sizeof(double *));
  for (i=0; i<numst; i++) {
    md->a[i]=(double *)calloc(numst,sizeof(double));
    if (fread(md->a[i], sizeof(double),numst, infile)!=numst)
      {
	fprintf(stderr, "**ERROR reading in model\n");
	return(0);
      }    
  }

  md->stpdf=(GaussModel **)calloc(numst, sizeof(GaussModel *));
  for (i=0; i<numst; i++)
    md->stpdf[i]=(GaussModel *)calloc(1, sizeof(GaussModel));

  for (i=0; i<numst; i++) {
    if (fread(&(md->stpdf[i]->exist), sizeof(int),1,infile)!=1) {
      fprintf(stderr, "**ERROR reading in model\n");
      return(0);
    }

    if (fread(&(md->stpdf[i]->dim), sizeof(int), 1,infile)!=1) {
      fprintf(stderr, "**ERROR reading in model\n");
      return(0);
    } 
    
    md->stpdf[i]->mean=(double *)calloc(md->stpdf[i]->dim,sizeof(double));
    if (fread(md->stpdf[i]->mean, sizeof(double), md->stpdf[i]->dim, 
	       infile)!=md->stpdf[i]->dim) {
      fprintf(stderr, "**ERROR reading in model\n");
      return(0);
    }
 
    if (fread(&(md->stpdf[i]->sigma_det), sizeof(double),1,infile)!=1) {
      fprintf(stderr, "**ERROR reading in model\n");
      return(0);
    }
 
    md->stpdf[i]->sigma=(double **)calloc(md->stpdf[i]->dim,sizeof(double *));
    for (m=0; m<md->stpdf[i]->dim; m++) {
      md->stpdf[i]->sigma[m]=(double *)calloc(md->stpdf[i]->dim, 
					      sizeof(double));
      if (fread(md->stpdf[i]->sigma[m], sizeof(double), md->stpdf[i]->dim,
		infile)!=md->stpdf[i]->dim) {
	fprintf(stderr, "**ERROR reading in model\n");
	return(0);
      } 
    }

    md->stpdf[i]->sigma_inv=(double **)calloc(md->stpdf[i]->dim,
					      sizeof(double *));
    for (m=0; m<md->stpdf[i]->dim; m++) {
      md->stpdf[i]->sigma_inv[m]=(double *)calloc(md->stpdf[i]->dim, 
						  sizeof(double));
      if (fread(md->stpdf[i]->sigma_inv[m], sizeof(double), md->stpdf[i]->dim,
		infile)!=md->stpdf[i]->dim) {
	fprintf(stderr, "**ERROR reading in model\n");
	return(0);
      } 
    }
  }
  
  return(1);
}


/*------------------------------------------------------------*/
unsigned char print_model(HmmModel *md, FILE *outfile)
{
  int i,j,k,m,n;
  int dim,numst,numcls;

  dim=md->dim;
  numst=md->numst;
  numcls=md->numcls;
  
  fprintf(outfile, "dim=%d\n", dim);
  fprintf(outfile, "numst=%d\n", numst);
  fprintf(outfile, "numcls=%d\n", numcls);

  fprintf(outfile, "\nState class membership:\n");
  for (i=0; i<numst; i++)
    fprintf(outfile, "%d ", md->stcls[i]);
  fprintf(outfile, "\n");

  fprintf(outfile, "\nTransition probability a00:\n");
  for (i=0; i<numst; i++)
    fprintf(outfile, "%8e ", md->a00[i]);
  fprintf(outfile, "\n");

  fprintf(outfile, "Transition probability a:\n");
  for (i=0; i<numst; i++) {
    for (j=0; j<numst; j++)
      fprintf(outfile, "%8e ", md->a[i][j]);
    fprintf(outfile, "\n");
  }

  fprintf(outfile, "\nThe Gaussian distributions of states:\n");
  for (i=0; i<numst; i++) {
    fprintf(outfile, "\nState %d =============================\n", i);
    fprintf(outfile, "exist=%d, dim=%d\n", md->stpdf[i]->exist, 
	    md->stpdf[i]->dim);

    fprintf(outfile, "Mean vector:\n");
    for (j=0; j<dim; j++)
      fprintf(outfile, "%.5e ", md->stpdf[i]->mean[j]);
    fprintf(outfile, "\n");

    fprintf(outfile, "Sigma_det=%e\n",md->stpdf[i]->sigma_det);

    fprintf(outfile, "Covariance matrix Sigma:\n");
 
    for (m=0; m<md->stpdf[i]->dim; m++) {
      for (n=0; n<md->stpdf[i]->dim; n++)
	fprintf(outfile, "%.5e ", md->stpdf[i]->sigma[m][n]);
      fprintf(outfile, "\n");
    }

    fprintf(outfile, "Covariance matrix inverse Sigma_inv:\n");
 
    for (m=0; m<md->stpdf[i]->dim; m++) {
      for (n=0; n<md->stpdf[i]->dim; n++)
	fprintf(outfile, "%.5e ", md->stpdf[i]->sigma_inv[m][n]);
      fprintf(outfile, "\n");
    }
  }
  
  return(1);
}

