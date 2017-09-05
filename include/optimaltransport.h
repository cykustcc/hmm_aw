//
//  optimaltransport.h
//  hmm_aw
//
//  Created by Yukun Chen on 4/27/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#ifndef optimaltransport_h
#define optimaltransport_h

double match_by_distmat_BADMM(int n,
                              int m,
                              double *C,
                              double* wX,
                              double* wY,
                              /** OUT **/ double *x,
                              /** OUT **/ double *lambda);

#endif /* optimaltransport_h */
