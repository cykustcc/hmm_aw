//
//  dist_utils.hpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/3/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#ifndef dist_utils_hpp
#define dist_utils_hpp

#include <stdio.h>
#include "hmm.h"

void pdist2(int d, size_t n, size_t m, double* A, double* B, double *C);

void calc_distmat(HmmModel* hmm1, HmmModel* hmm2);
#endif /* dist_utils_hpp */
