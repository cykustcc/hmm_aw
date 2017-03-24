//
//  mosek_solver.h
//  hmm_aw
//
//  Created by Yukun Chen on 3/8/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#ifndef mosek_solver_h
#define mosek_solver_h
#include <vector>
#ifdef __cplucplus
extern "C"{
#endif
  void solver_setup();
  void solver_release();
  double match_by_distmat(int n,
                          int m,
                          double *C,
                          std::vector<double> &wX,
                          std::vector<double> &wY,
                          /** OUT **/ double *x,
                          /** OUT **/ double *lambda);
#ifdef __cplucplus
}
#endif


#endif /* mosek_solver_h */
