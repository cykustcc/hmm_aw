#ifndef __hmm_aw_cluster_h
#define __hmm_aw_cluster_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern void split(float *cdwd, float *newcdwd, int dim, float *stdv);
extern void centroid(float *cdbk, int dim, int numcdwd, float *vc,
		int *index, int numdata);
extern void cellstdv(float *cdbk, float *stdv, int dim, int numcdwd, 
		     float *vc, int *index,  int numdata);
extern float mse_dist(float *cdwd, float *vc, int dim);
extern void encode(float *cdbk, int dim, int numcdwd, float *vc, int *code,
		   int numdata);
extern float lloyd(float *cdbk, int dim, int numcdwd, float *vc, int numdata, 
		   float threshold);
extern float kmeans(float *cdbk, int dim, int maxnumcdwd, int *fnumcdwd, 
		    float *vc, int numdata, float threshold, float distthred);

#endif

