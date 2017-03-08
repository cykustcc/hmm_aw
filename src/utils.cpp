#include "utils.h"
#include <stdio.h>

void matInit(float* X, unsigned int size, float value)
{
  for ( int i = 0 ; i < size ; i++  ){
    X[i] = value;
  }
}

void print_mat(float* X, int dim, int n){
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("%f,\t", X[j+i*dim]);
    }
      printf(";\n");
  }
}
