#include "utils.h"
#include <stdio.h>
#include <iostream>

void matInit(float* X, unsigned int size, float value)
{
  for ( int i = 0 ; i < size ; i++  ){
    X[i] = value;
  }
}

void print_mat_float(float* X, int dim, int n){
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("%f,\t", X[j+i*dim]);
    }
      printf(";\n");
  }
}

void print_mat_float(std::vector<float>& X, int dim, int n){
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("%f,\t", X[j+i*dim]);
    }
    printf(";\n");
  }
}

void print_mat_double(double* X, int dim, int n){
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("%lf,\t", X[j+i*dim]);
    }
    printf(";\n");
  }
}

void print_2d_double(std::vector<std::vector<double>>& a){
  size_t n = a.size();
  size_t m = n > 0? a[0].size(): 0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      std::cout<<a[i][j]<<" ";
    }
    std::cout<<std::endl;
  }
}
