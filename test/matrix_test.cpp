//
//  matrix_test.cpp
//  hmm_aw
//
//  Created by Yukun Chen on 3/1/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include <iostream>
#include "matrix.h"
#include "gtest/gtest.h"

using namespace std;
class MatrixTest: public ::testing::Test{
protected:
  std::vector<std::vector<float>> P = {{4,12,-16},
                                      {12,37,-43},
                                      {-16,-43,98}};
};

TEST_F(MatrixTest, CholeskyDecomp){
  std::vector<std::vector<float>> S(3, std::vector<float>(3, 0.0));
  cholesky_decomp(P, S);
  print_matrix(S);
}
