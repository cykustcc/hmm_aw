#include "utils.h"
#include "matrix.h"
#include "gtest/gtest.h"
#include <glog/logging.h>
#include <gflags/gflags.h>

class UtilsTest: public ::testing::Test{
protected:
    float mat_1x3[3] = {1.0, 1.0, 1.0};
    float mat_1x10[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float mat_5x5[5][5] = {
        {1.0 , 0.0, 0.0, 0.0, 0.0},
        {0.0 , 1.0, 0.0, 0.0, 0.0},
        {0.0 , 0.0, 1.0, 0.0, 0.0},
        {0.0 , 0.0, 0.0, 1.0, 0.0},
        {0.0 , 0.0, 0.0, 0.0, 1.0}
        };
};
