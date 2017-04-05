//
//  main.cpp
//  testlapacke
//
//  Created by Yukun Chen on 4/3/17.
//  Copyright © 2017 cyk. All rights reserved.
//

#include <iostream>
#include <cassert>
#include "lapacke.h"
#include <vector>

#include <stdio.h>
// #include "lapacke_example_aux.h"
void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
  lapack_int i, j;
  printf( "\n %s\n", desc );
  
  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i*ldm+j] );
    printf( "\n" );
  }
}


/* Auxiliary routine: printing a matrix */
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
  lapack_int i, j;
  printf( "\n %s\n", desc );
  
  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i+j*ldm] );
    printf( "\n" );
  }
}

/* Auxiliary routine: printing a vector of integers */
void print_vector( char* desc, lapack_int n, lapack_int* vec ) {
  lapack_int j;
  printf( "\n %s\n", desc );
  for( j = 0; j < n; j++ ) printf( " %6i", vec[j] );
  printf( "\n" );
}


/* Main program */
void example1()
{
  /* Locals */
//  double A[5][3] = {1,1,1,2,3,4,3,5,2,4,2,5,5,4,3};
//  double b[5][2] = {-10,-3,12,14,14,12,16,16,18,16};
//  std::vector<double> A({1,1,1,2,3,4,3,5,2,4,2,5,5,4,3});
//  std::vector<double> b({-10,-3,12,14,14,12,16,16,18,16});
  std::vector<std::vector<double>> A ={{1,1,1,2,3,4,3,5,2,4,2,5,5,4,3}};
  std::vector<std::vector<double>> b ={{-10,-3,12,14,14,12,16,16,18,16}};
//  std::vector<std::vector<double>> A = {{1,1,1},{2,3,4},{3,5,2},{4,2,5},{5,4,3}};
//  std::vector<std::vector<double>> b = {{-10,-3},{12,14},{14,12},{16,16},{18,16}};
  lapack_int info,m,n,lda,ldb,nrhs;
  int i,j;
  
  /* Initialization */
  m = 5;
  n = 3;
  nrhs = 2;
  lda = 3;
  ldb = 2;
  
  /* Print Entry Matrix */
//  print_matrix_rowmajor( "Entry Matrix A", m, n, *A, lda );
  print_matrix_rowmajor( "Entry Matrix A", m, n, &A[0][0], lda );
  /* Print Right Rand Side */
//  print_matrix_rowmajor( "Right Hand Side b", n, nrhs, *b, ldb );
  print_matrix_rowmajor( "Right Hand Side b", n, nrhs, &b[0][0], ldb );
  printf( "\n" );
  
  /* Executable statements */
  printf( "LAPACKE_dgels (row-major, high-level) Example Program Results\n" );
  /* Solve least squares problem*/
//  info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,*A,lda,*b,ldb);
  info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,&A[0][0],lda,&b[0][0],ldb);
  
  /* Print Solution */
//  print_matrix_rowmajor( "Solution", n, nrhs, *b, ldb );
  print_matrix_rowmajor( "Solution", n, nrhs, &b[0][0], ldb );
  printf( "\n" );
  exit( 0 );
} /* End of LAPACKE_dgels Example */


double A[36] = {4, 3, 0, 4, 0, 3, 2, 3, 3, 2, 2, 0, 0, 3, 2, 2, 0, 3,
  3, 3, 2, 3, 0, 0, 1, 2, 3, 4, 0, 2, 2, 4, 2, 4, 3, 2};

double Arow[36] = {
  4, 2, 0, 3, 1, 2, 3, 3, 3, 3, 2, 4, 0, 3, 2, 2, 3, 2,
  4, 2, 2, 3, 4, 4, 0, 2, 0, 0, 0, 3, 3, 0, 3, 0, 2, 2,
};

double B[24] = {2, 2, 3, 2, 4, 0, 3, 1, 3, 2, 1, 0,
  2, 4, 3, 4, 2, 0, 1, 4, 2, 4, 0, 4};

double Brow[24] = {
  2, 3, 2, 1, 2, 1, 4, 4, 3, 3, 3, 2, 2, 2, 4, 4, 4, 1, 2, 0, 0, 0, 0, 4,
};

/* Auxiliary routine: printing a matrix */
void print_matrix(char *desc, char major, int m, int n, double *a, int lda) {
  int i, j;
  printf("\n %s (%c major)\n", desc, major);
  if (major == 'C')
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++)
        printf(" %6.2f", a[i + j * lda]);
      printf("\n");
    }
  
  else {
    printf("m=%d,n=%d\n", m, n);
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++)
        printf(" %6.2f", a[j + i * lda]);
      printf("\n");
    }
  }
}

void eq_solve(double *a, double *b, int n, int nrhs, char major) {
  int info, i;
  // char trans = 'n';
  lapack_int lda, ldb, *ipiv;
  lda = n;
  ipiv = (int *)malloc(sizeof(int) * n);
  /* Print A */
  print_matrix("A", major, n, n, a, lda);
  
  if (major == 'C') {
    ldb = n;
    /* Print B */
    print_matrix("B", major, n, nrhs, b, ldb);
    // EQUIVALENT STATEMENT
    info = LAPACKE_dgesv(LAPACK_COL_MAJOR, (lapack_int)n, (lapack_int)nrhs, a,
                         lda, ipiv, b, ldb);
    // dgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
  } else {
    ldb = nrhs;
    /* Print B */
    print_matrix("B", major, n, nrhs, b, ldb);
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int)n, (lapack_int)nrhs, a,
                         lda, ipiv, b, ldb);
  }
  
  if (info != 0) {
    fprintf(stderr, "dgesv: info = %d\n", info);
  }
  assert(info == 0);
  printf("dgesv passed\n");
  // debug_print("b using lapacke", b,n,nrhs,0);
  // print answer here
  print_matrix("Answer", major, n, nrhs, b, ldb);
  
  free(ipiv);
}

int main(void) {
//  printf("==============================\n");
//  printf("\tSOLVING IN COL MAJOR\n");
//  printf("==============================\n");
//  eq_solve(A, B, 6, 4, 'C');
//  printf("\n==============================\n");
//  printf("\tSOLVING IN ROW MAJOR\n");
//  printf("==============================\n");
//  eq_solve(Arow, Brow, 6, 4, 'R');
  example1();
  return 0;
}
