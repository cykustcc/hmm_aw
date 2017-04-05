//
//  test_common.h
//  hmm_aw
//
//  Created by Yukun Chen on 3/11/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#ifndef test_common_h
#define test_common_h
#include<string>
#include "hmm.h"


extern std::string root_path;


void expect_same_hmm(HmmModel& md1, HmmModel& md2, double tol);

#endif /* test_common_h */
