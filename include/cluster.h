#ifndef __hmm_aw_cluster_h
#define __hmm_aw_cluster_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

extern void split(std::vector<float> &cdwd,
                  int cdwd_baseidx,
                  int newcdwd_baseidx,
                  int dim,
                  std::vector<float> &stddev,
                  int stddev_baseidx);

extern void centroid(std::vector<float> &cdbk,
                     int dim,
                     int numcdwd,
                     std::vector<float> &vc,
                     std::vector<int> &index,
                     int numdata);

extern void cellstdv(std::vector<float> &cdbk,
                     std::vector<float> &stdv,
                     int dim,
                     int numcdwd,
                     std::vector<float> &vc,
                     std::vector<int> &index,
                     int numdata);

extern float mse_dist(std::vector<float> &cdwd,
                      std::vector<float> &vc,
                      int dim);

extern void encode(std::vector<float> &cdbk,
                   int dim,
                   int numcdwd,
                   std::vector<float> &vc,
                   std::vector<int> &code,
                   int numdata);

extern float lloyd(std::vector<float> &cdbk,
                   int dim,
                   int numcdwd,
                   std::vector<float> &vc,
                   int numdata,
                   float threshold);

extern float kmeans(std::vector<float> &cdbk,
                    int dim, int maxnumcdwd,
                    float &fnumcdwd,
                    std::vector<float> &vc,
                    int numdata,
                    float threshold,
                    float distthred);

#endif

