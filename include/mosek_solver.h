//
//  mosek_solver.h
//  hmm_aw
//
//  Created by Yukun Chen on 3/8/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#ifndef mosek_solver_h
#define mosek_solver_h

#ifdef __cplucplus
extern "C"{
#endif
  double match_by_distmat(int n, int m, double *C, double *wX, double *wY,
                             /** OUT **/ double *x, /** OUT **/ double *lambda);

#ifdef __cplucplus
}
#endif


#endif /* mosek_solver_h */
