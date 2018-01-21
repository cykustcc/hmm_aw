#ifndef __hmm_aw_matrix_h
#define __hmm_aw_matrix_h
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
/*-------------------------- mat_simple.c --------------------------*/

/*-------------------------------------------------------------------*/
/*------------------------ Print out Vector -------------------------*/
/*-------------------------------------------------------------------*/

template<typename DType>
void print_vector(std::vector<DType>& vt) {
  int size = vt.size();
  for (int i = 0; i < size; i++) {
    std::cout<<vt[i]<<" ";
  }
  std::cout<<std::endl;
}

template<typename DType>
void print_vector(std::vector<DType>& vt, int m, int n) {
  int size = vt.size();
  assert(m * n <= size);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << vt[i * n + j] << " ";
    }
    std::cout << std::endl;
  }
}


/*-------------------------------------------------------------------*/
/*------------------------ Print out Matrix -------------------------*/
/*-------------------------------------------------------------------*/
template<typename DType>
extern void print_matrix(std::vector<std::vector<DType>> &mt) {
  int rows = (int) mt.size();
  int cols = rows > 0 ? (int) mt[0].size(): 0;

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      std::cout << mt[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

template<typename DType>
void print_matrix(DType* X, int dim, int n) {
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("%lf,\t", X[j + i * dim]);
    }
    printf(";\n");
  }
}

/*-------------------------------------------------------------------*/
/*---------------      Numerical Programs             ---------------*/
/*---------------      LU decomposition programs      ---------------*/
/*---------------      Calculate matrix inversion     ---------------*/
/*---------------      Calculate matrix determinant   ---------------*/
/*-------------------------------------------------------------------*/
template<typename DType>
extern DType mat_det(std::vector<std::vector<DType>> mt,
                     int dim) {
  DType res;

  if (dim == 1)
    return(mt[0][0]);

  std::vector<std::vector<DType>> submt(dim - 1,
                                        std::vector<DType>(dim - 1, 0.0));

  for (int i = 1; i < dim; i++) {
    for (int j = 1; j < dim; j++) {
      submt[i - 1][j - 1] = mt[i][j];
    }
  }

  int n = 1;
  res = 0.0;

  for (int i = 0; i < dim; i++) {
    res += (n*mt[i][0] * mat_det(submt, dim - 1));
    n = -n;
    if (i == dim - 1)
      continue;
    for (int j = 1; j < dim; j++) {
      submt[i][j-1] = mt[i][j];
    }
  }

  return(res);
}

template<typename DType>
unsigned char ludcmp(std::vector<std::vector<DType>> &a,
                     int n,
                     std::vector<int> &indx,
                     DType &d)
{
  int imax = 0;
  DType big, dum, sum, temp;
  std::vector<float> vv(n, 0.0);
  DType TINY = 1e-20;

  d = 1.0;
  for (int i = 0; i < n; i++) {
    big = 0.0;
    for (int j = 0; j < n; j++)
      if ((temp = fabs(a[i][j])) > big)
        big = temp;
    if (big == 0.0) {
      fprintf(stderr, "Singular matrix in ludcmp\n");
      return(2); /* 2 stands for singular matrix */
    }
    vv[i]=1.0 / big;
  }

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < j; i++) {
      sum = a[i][j];
      for (int k = 0; k < i; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (int i = j; i < n; i++) {
      sum = a[i][j];
      for (int k = 0; k < j; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum=vv[i] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (int k = 0; k < n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0) a[j][j] = TINY;
    if (j != n-1) {
      dum = 1.0 / (a[j][j]);
      for (int i = j+1; i < n; i++)
        a[i][j] *= dum;
    }
  }
  return(1);
}

//	Cholesky_Decomposition returns the Cholesky Decomposition Matrix.
// mt_P = mt_S * mt_S^T
template<typename DType>
extern void cholesky_decomp(std::vector<std::vector<DType>> &mt_P,
                            /*Output: */ std::vector<std::vector<DType>> &mt_S,
                            bool diag=false)
{
  DType temp = 0, temp2 = 0;
  size_t m = mt_P.size(), n = m>0?mt_P[0].size():0;
  assert(m == n && mt_S.size() == m);
  //	Initialize and populate matrix L which will be the lower Cholesky
  if (diag) {
    for (int i = 0; i < m; i++){
      mt_S[i][i] = sqrt(mt_P[i][i]);
    }
  }else{
    for (int i = 0; i < m; i++){
      for (int j = 0; j < m; j++){
        temp = 0; temp2 = 0;
        if (i > j){
          if (j > 0){
            for (int k = 1; k < j + 1; k++)
              temp2 += (mt_S[i][k - 1] * mt_S[j][k - 1]);
          }
          mt_S[i][j] = (mt_P[i][j] - temp2) / mt_S[j][j];
        }
        else if (i == j){
          for (int k = 0; k < i; k++)
            temp += pow(mt_S[i][k], 2);
          mt_S[i][j] = sqrt(mt_P[i][j] - temp);
        }
        else
          mt_S[i][j] = 0;
      }
    }
  }
}

template<typename DType>
extern void lubksb(std::vector<std::vector<DType>> &a,
                   int n,
                   std::vector<int> &indx,
                   std::vector<DType> &b){
  int ii = -1, ip;
  DType sum;

  for (int i = 0; i<n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii >= 0)
      for (int j = ii; j < i; j++) sum -= a[i][j]*b[j];
    else if (sum != 0.0) ii = i;
    b[i] = sum;
  }
  for (int i = n-1; i >= 0; i--) {
    sum = b[i];
    for (int j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

template<typename DType>
unsigned char mat_inv(std::vector<std::vector<DType>> &mt,
                      std::vector<std::vector<DType>> &y,
                      int dim){
  DType d;
  std::vector<std::vector<DType>> a(dim, std::vector<DType>(dim, 0.0));

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      a[i][j] = mt[i][j];
    }
  }

  std::vector<DType> col(dim);
  std::vector<int> indx(dim);

  ludcmp(a, dim, indx, d);
  for (int j = 0; j < dim; j++) {
    for (int i = 0; i < dim; i++) col[i] = 0.0;
    col[j] = 1.0;
    lubksb(a, dim, indx, col);
    for (int i = 0; i < dim;i++)
      y[i][j] = col[i];
  }

  return(1);
}

template<typename DType>
float mat_det_ludcmp(std::vector<std::vector<DType>> &mt,
                     int dim)
{
  DType d;

  std::vector<std::vector<DType>> a(dim, std::vector<DType>(dim, 0.0));
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      a[i][j] = mt[i][j];
    }
  }
  std::vector<int> indx(dim, 0);

  ludcmp(a, dim, indx, d);
  for (int j = 0; j < dim; j++) d *= a[j][j];

  return(d);
}




/** compute the determinant and at the same provide the inverse matrix **/
/** combining these two operations save computation since they share   **/
/** the call to ludcmp_double                      **/
template<typename DType>
unsigned char mat_det_inv(std::vector<std::vector<DType>> &mt,
                          std::vector<std::vector<DType>> &y,
                          DType &det,
                          int dim)
{
  double d;
  std::vector<std::vector<DType>> a(mt);

  /** initialize matrix determinant **/
  det = 0.0;
  std::vector<DType> col(dim, 0.0);
  std::vector<int> indx(dim, 0);

  int m = ludcmp(a, dim, indx, d);
  if (m == 2) {  /** singular matrix **/
    det = 0.0;
    return(2);
  }

  for (int j = 0; j < dim; j++) d *= a[j][j];
  det = d;

  for (int j = 0; j < dim; j++) {
    for (int i = 0; i < dim; i++) col[i] = 0.0;
    col[j] = 1.0;
    lubksb(a, dim, indx, col);
    for (int i = 0; i < dim; i++)
      y[i][j] = col[i];
  }
  return(1);
}

#endif  //__hmm_aw_matrix_h
