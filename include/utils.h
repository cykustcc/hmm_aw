#ifndef UTILS_H_HMM_AW
#define UTILS_H_HMM_AW
#include <vector>

void matInit(float* X, unsigned int size, float value);
void print_mat_float(float* X, int dim, int n);
void print_mat_float(std::vector<float>& X, int dim, int n);
void print_mat_double(double* X, int dim, int n);

void print_2d_double(std::vector<std::vector<double>>& a);

#endif
