#include "cluster.h"

void split(float *cdwd, float *newcdwd, int dim, float *stddev)
{
  float mult_offset=0.1;
  int i;

  /** if cdwd is way off zero and the variance is small     **/
  /** the offset may be out of the data range and generates **/
  /** empty cells.                                          **/
  /*****
  for (i=0; i<dim; i++) {
    newcdwd[i] = cdwd[i]*(1+mult_offset*drand48());
  }
  *****/

  /* set random range to [0.25, 0.75] */
  for (i=0; i<dim; i++) {
    newcdwd[i] = cdwd[i]+stddev[i]*mult_offset*(0.25+drand48()/2.0);
  }
    
}

void centroid(float *cdbk, int dim, int numcdwd, float *vc,
	 int *index, int numdata)
{
  int i,j,k,m,n;
  int *ct;

  ct=(int *)malloc(numcdwd*sizeof(int));

  if (index==NULL)
    {
      for (k=0; k<dim; k++)
	cdbk[k] = 0.0;
      for (i=0; i<numdata; i++) 
	for (k=0; k<dim; k++) 
	  cdbk[k]+=vc[i*dim+k];
      for (k=0; k<dim; k++)
	cdbk[k] /= ((float)numdata);
    }
  else
    {
      for (j=0; j<numcdwd; j++) {
	for (k=0; k<dim; k++)
	  cdbk[j*dim+k] = 0.0;
	ct[j] = 0;
      }
      
      for (i=0; i<numdata; i++) {
	for (k=0; k<dim; k++) 
	  cdbk[index[i]*dim+k]+=vc[i*dim+k];
	ct[index[i]]++;
      }
      
      for (j=0; j<numcdwd; j++) 
	for (k=0; k<dim; k++)
	  cdbk[j*dim+k] /= ((float)ct[j]);
    }

  free(ct);

}  

void cellstdv(float *cdbk, float *stddev, int dim, int numcdwd, float *vc,
	      int *index,  int numdata)
{
  int i,j,k,m,n;
  int *ct;

  ct=(int *)malloc(numcdwd*sizeof(int));

  for (j=0; j<numcdwd; j++) {
    for (k=0; k<dim; k++)
      stddev[j*dim+k] = 0.0;
    ct[j] = 0;
  }
      
  for (i=0; i<numdata; i++) {
    for (k=0; k<dim; k++) 
      stddev[index[i]*dim+k]+=((vc[i*dim+k]-cdbk[index[i]*dim+k])*
	(vc[i*dim+k]-cdbk[index[i]*dim+k]));
    ct[index[i]]++;
  }
  
  for (j=0; j<numcdwd; j++) { 
    if (ct[j]>0) {
      for (k=0; k<dim; k++) {
	stddev[j*dim+k] /= ((float)ct[j]);
	stddev[j*dim+k]=sqrt(stddev[j*dim+k]);
      }
    }
    else {
      for (k=0; k<dim; k++) stddev[j*dim+k]=1.0;
    }
  }

  free(ct);

}  


float mse_dist(float *cdwd, float *vc, int dim)
{
  float mse= 0.0;
  int i;

  for (i=0; i<dim; i++)
    mse += (vc[i]-cdwd[i])*(vc[i]-cdwd[i]);

  return(mse);
}


void encode(float *cdbk, int dim, int numcdwd, float *vc, int *code,
	    int numdata)
{
  int i,j,k,m,n;
  float *dist,minv;

  dist=(float *)calloc(numcdwd,sizeof(float));

  for (i=0; i<numdata; i++) {
    for (j=0; j<numcdwd;j++)
      dist[j]=mse_dist(cdbk+j*dim, vc+i*dim, dim);    
    code[i]=0;
    minv=dist[0];
    for (j=1;j<numcdwd;j++)
      if (dist[j]<minv) {
	minv=dist[j];
	code[i]=j;
      }
  }

  free(dist);
}


float lloyd(float *cdbk, int dim, int numcdwd, float *vc, int numdata, 
	    float threshold)
     /* return the value of the mean squared distance, i.e., average */
     /* squared distance between a sample and its centroid           */
     /* threshold is for controling when to stop the loop */
{
  int i,j,k,m,n;
  int ite;
  float dist, olddist, minmse, mse;
  int min_iteration=2;
  /*float threshold = 0.005;*/
  int *index;
  float *tp;
  int cdbksz2, temp_int, new_cdwds, cdbksz;
  float *stddev; // standard deviation for appropriate split

  cdbksz2 = 0;
  temp_int = 1;
  while (temp_int < numcdwd) {
    cdbksz2++;
    temp_int += temp_int;
  }  


  index = (int *)calloc(numdata, sizeof(int));
  stddev = (float *)calloc(numcdwd*dim,sizeof(float));

  centroid(cdbk, dim, 1, vc, NULL, numdata);

  /* compute standard deviation for each cell */
  for (i=0;i<numdata;i++) index[i]=0;
  cellstdv(cdbk,stddev,dim,numcdwd,vc,index,numdata);

  if (numcdwd==1) {
    dist = 0.0;
    for (i=0, k=0; i<numdata; i++)
      for (j=0; j<dim; j++, k++)
	dist += (cdbk[j]-vc[k])*(cdbk[j]-vc[k]);
    dist /= ((float)numdata);
  }

  cdbksz = 1;

  for (ite=0; ite<cdbksz2; ite++) {
    new_cdwds = ((numcdwd - 2*cdbksz) >= 0) ? cdbksz : numcdwd - cdbksz;

    for (k=0; k<new_cdwds; k++)
      split(cdbk+k*dim, cdbk+(cdbksz+k)*dim, dim, stddev+k*dim);

    cdbksz += new_cdwds;

    dist=HUGE;
    m = 0;

    while (m < min_iteration || 
	   (fabs((double)(olddist-dist))/olddist > threshold
	    && dist < olddist))
      {
	m++;
	olddist = dist;
	tp = vc;
	dist = 0.0;
	for (i=0; i<numdata; i++)
	  {
	    minmse = mse_dist(cdbk, tp, dim);
	    index[i]= 0;
	    
	    for (j=1; j<cdbksz; j++)
	      {
		mse = mse_dist(cdbk+j*dim, tp, dim);
		if (mse<minmse)
		  {
		    minmse=mse;
		    index[i]=j;
		  }
	      }
	    
	    dist += minmse;
	    tp += dim;
	  }
	
	dist /= ((float)numdata);

	centroid(cdbk, dim, cdbksz, vc, index, numdata);
      }
    cellstdv(cdbk,stddev,dim,cdbksz,vc,index,numdata);

  }

  free(index);
  free(stddev);

  return(dist);

}	    


/** The kmeans algorithm is close to lloyd, except that the number **/
/** codewords is not preselected.  Instead, it is determined by    **/
/** the minimum number of codewords that lead to a mean squared    **/
/** error below a given threshold.  The number of codewords is     **/
/** upper bounded by the given maxnumcdwd.                         **/

float kmeans(float *cdbk, int dim, int maxnumcdwd, int *fnumcdwd, 
	     float *vc, int numdata, float threshold, float distthred)
{
  int i,j,k,m,n;
  int ite, splitwd;
  float dist, olddist, minmse, mse;
  int min_iteration=2;
  int numcdwd;

  centroid(cdbk, dim, numcdwd, vc, NULL, numdata);

  dist = 0.0;
  for (i=0, k=0; i<numdata; i++)
    for (j=0; j<dim; j++, k++)
      dist += (cdbk[j]-vc[k])*(cdbk[j]-vc[k]);
  dist /= ((float)numdata);
  if (dist<distthred) {
    *fnumcdwd=1;
    return(dist);
  }

  numcdwd=2;
  do {
    dist=lloyd(cdbk, dim, numcdwd, vc, numdata, threshold);
    numcdwd++;
    //fprintf(stderr, "numcdwd=%d, dist=%f\n", numcdwd,dist);
  } while (numcdwd<=maxnumcdwd && dist > distthred);
  
  *fnumcdwd=numcdwd-1;
  
  return(dist);
}	    

