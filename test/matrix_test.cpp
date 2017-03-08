//
//  matrix_test.cpp
//  hmmjia
//
//  Created by Yukun Chen on 3/1/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include <iostream>
#include "matrix.h"
#include "gtest/gtest.h"

using namespace std;
class VectorTest: public ::testing::Test{
protected:
    float mat_1x3[3] = {1.0, 1.0, 1.0};
    float mat_1x10[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
};

TEST(GeneralVector, TestVectorUchar){
    cout<<"hello"<<endl;
}
