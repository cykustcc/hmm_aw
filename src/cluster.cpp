#include "cluster.h"
#include "glog/logging.h"

void split(std::vector<float> &cdwd,
           int cdwd_baseidx,
           int newcdwd_baseidx,
           int dim,
           std::vector<float> &stddev,
           int stddev_baseidx)
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
    cdwd[newcdwd_baseidx + i] = cdwd[cdwd_baseidx + i]+stddev[stddev_baseidx + i]*mult_offset*(0.25+drand48()/2.0);
  }
}

void centroid(std::vector<float> &cdbk,
              int dim, int numcdwd,
              std::vector<float> &vc,
              std::vector<int> &index,
              int numdata)
{
  std::vector<int> ct(numcdwd,0);

  if (index.size()==0){
    for (int k=0; k<dim; k++)
      cdbk[k] = 0.0;
    for (int i=0; i<numdata; i++)
      for (int k=0; k<dim; k++)
        cdbk[k]+=vc[i*dim+k];
    for (int k=0; k<dim; k++)
      cdbk[k] /= ((float)numdata);
  }else{
    for (int j=0; j<numcdwd; j++) {
      for (int k=0; k<dim; k++)
        cdbk[j*dim+k] = 0.0;
      ct[j] = 0;
    }
      
    for (int i=0; i<numdata; i++) {
      for (int k=0; k<dim; k++)
        cdbk[index[i]*dim+k]+=vc[i*dim+k];
      ct[index[i]]++;
    }
      
    for (int j=0; j<numcdwd; j++)
      for (int k=0; k<dim; k++)
        cdbk[j*dim+k] /= ((float)ct[j]);
    }
}  

void cellstdv(std::vector<float> &cdbk,
              std::vector<float> &stddev,
              int dim,
              int numcdwd,
              std::vector<float> &vc,
              std::vector<int> &index,
              int numdata)
{
  std::vector<int> ct(numcdwd, 0);

  for (int j=0; j<numcdwd; j++) {
    for (int k=0; k<dim; k++)
      stddev[j*dim+k] = 0.0;
    ct[j] = 0;
  }
      
  for (int i=0; i<numdata; i++) {
    for (int k=0; k<dim; k++)
      stddev[index[i]*dim+k]+=((vc[i*dim+k]-cdbk[index[i]*dim+k])*
      (vc[i*dim+k]-cdbk[index[i]*dim+k]));
    ct[index[i]]++;
  }
  
  for (int j=0; j<numcdwd; j++) {
    if (ct[j]>0) {
      for (int k=0; k<dim; k++) {
	stddev[j*dim+k] /= ((float)ct[j]);
	stddev[j*dim+k]=sqrt(stddev[j*dim+k]);
      }
    }
    else {
      for (int k=0; k<dim; k++) stddev[j*dim+k]=1.0;
    }
  }
}  


float mse_dist(std::vector<float> &cdwd,
               int cdwd_baseidx,
               std::vector<float> &vc,
               int vc_baseidx,
               int dim)
{
  float mse= 0.0;

  for (int i=0; i<dim; i++)
    mse += (vc[vc_baseidx+i]-cdwd[cdwd_baseidx+i])*(vc[vc_baseidx+i]-cdwd[cdwd_baseidx+i]);

  return(mse);
}


void encode(std::vector<float> &cdbk,
            int dim,
            int numcdwd,
            std::vector<float> &vc,
            std::vector<int> &code,
            int numdata){
  float minv;
  std::vector<float> dist(numcdwd,0.0);

  for (int i=0; i<numdata; i++) {
    for (int j=0; j<numcdwd;j++)
      dist[j]=mse_dist(cdbk, j*dim, vc, i*dim, dim);
    code[i]=0;
    minv=dist[0];
    for (int j=1;j<numcdwd;j++)
      if (dist[j]<minv) {
        minv=dist[j];
        code[i]=j;
      }
  }
}


float lloyd(std::vector<float> &cdbk,
            int dim,
            int numcdwd,
            std::vector<float> &vc,
            int numdata,
            float threshold)
     /* return the value of the mean squared distance, i.e., average */
     /* squared distance between a sample and its centroid           */
     /* threshold is for controling when to stop the loop */
{
  float dist = 0.0, olddist = 0.0, minmse, mse;
  int min_iteration=2;
  /*float threshold = 0.005;*/
  int cdbksz2, temp_int, new_cdwds, cdbksz;

  cdbksz2 = 0;
  temp_int = 1;
  while (temp_int < numcdwd) {
    cdbksz2++;
    temp_int += temp_int;
  }  

  std::vector<int> index(numdata, 0);
  std::vector<float> stddev(numcdwd*dim, 0.0);// standard deviation for appropriate split
  std::vector<int> tmpindex;
  
  centroid(cdbk, dim, 1, vc, tmpindex, numdata);

  /* compute standard deviation for each cell */
  for (int i=0;i<numdata;i++) index[i]=0;
  cellstdv(cdbk,stddev,dim,numcdwd,vc,index,numdata);

  if (numcdwd==1) {
    dist = 0.0;
    for (int i=0, k=0; i<numdata; i++)
      for (int j=0; j<dim; j++, k++)
          dist += (cdbk[j]-vc[k])*(cdbk[j]-vc[k]);
    dist /= ((float)numdata);
  }

  cdbksz = 1;

  for (int ite=0; ite<cdbksz2; ite++) {
    new_cdwds = ((numcdwd - 2*cdbksz) >= 0) ? cdbksz : numcdwd - cdbksz;

    for (int k=0; k<new_cdwds; k++)
      split(cdbk, k*dim, (cdbksz+k)*dim, dim, stddev, k*dim);

    cdbksz += new_cdwds;

    dist=HUGE;
    int m = 0;

    while (m < min_iteration || 
	    (fabs((double)(olddist-dist))/olddist > threshold
	    && dist < olddist)){
      m++;
      olddist = dist;
      dist = 0.0;
      for (int i=0; i<numdata; i++){
        minmse = mse_dist(cdbk, 0, vc, i*dim, dim);
        index[i]= 0;
        for (int j=1; j<cdbksz; j++){
          mse = mse_dist(cdbk, j*dim, vc, i*dim, dim);
          if (mse<minmse){
            minmse=mse;
            index[i]=j;
          }
        }
        dist += minmse;
      }
      dist /= ((float)numdata);

      centroid(cdbk, dim, cdbksz, vc, index, numdata);
    }
    cellstdv(cdbk,stddev,dim,cdbksz,vc,index,numdata);
  }

  return(dist);
}	    


/** The kmeans algorithm is close to lloyd, except that the number **/
/** codewords is not preselected.  Instead, it is determined by    **/
/** the minimum number of codewords that lead to a mean squared    **/
/** error below a given threshold.  The number of codewords is     **/
/** upper bounded by the given maxnumcdwd.                         **/

float kmeans(std::vector<float> &cdbk,
             int dim,
             int maxnumcdwd,
             float &fnumcdwd,
             std::vector<float> &vc,
             int numdata,
             float threshold,
             float distthred)
{
  int ite, splitwd;
  float dist, olddist, minmse, mse;
  int min_iteration=2;
  int numcdwd;

  std::vector<int> index;
  centroid(cdbk, dim, numcdwd, vc, index, numdata);

  dist = 0.0;
  
  for (int i=0, k=0; i<numdata; i++)
    for (int j=0; j<dim; j++, k++)
      dist += (cdbk[j]-vc[k])*(cdbk[j]-vc[k]);
  
  dist /= ((float)numdata);
  if (dist<distthred) {
    fnumcdwd=1;
    return(dist);
  }

  numcdwd=2;
  do {
    dist=lloyd(cdbk, dim, numcdwd, vc, numdata, threshold);
    numcdwd++;
    //fprintf(stderr, "numcdwd=%d, dist=%f\n", numcdwd,dist);
  } while (numcdwd<=maxnumcdwd && dist > distthred);
  
  fnumcdwd=numcdwd-1;
  
  return(dist);
}	    

