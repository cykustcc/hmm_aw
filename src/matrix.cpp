/**************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*                                                                        */
/*                        Jia  Li                                         */
/*                        April 2, 1998                                   */
/*                                                                        */
/*--------------------------------------------------------------------------
 ^^------------------------------------------------------------------------^^
 ^^                                                                        ^^
 ^^  matrix.c                                                              ^^
 ^^------------------------------------------------------------------------^^
 ^^                                                                        ^^
 ^^   Function:                                                            ^^
 ^^   Subroutines for dynamic memory operation, including assignment,      ^^
 ^^   free memory and copy.                                                ^^
 ^^                                                                        ^^
 ^^------------------------------------------------------------------------^^
 ^^                                                                        ^^
 ^^   Includes:                                                            ^^
 ^^   ......                                                               ^^
 ^^                                                                        ^^
 ^^------------------------------------------------------------------------^^
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/**************************************************************************/


#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include <stdlib.h>
#include "glog/logging.h"
#include <vector>




/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*------------------------ Print out Matrix -------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void print_matrix_uchar(std::vector<std::vector<unsigned char>> &mt,
                        int rows,
                        int cols)
{
  for (int i=0; i<rows; i++) {
    for (int j=0; j<cols; j++) {
      fprintf(stdout, "%d  ", mt[i][j]);
      if ((j+1)%8==0)
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }
}

void print_matrix_int(std::vector<std::vector<int>> &mt,
                      int rows,
                      int cols)
{
  for (int i=0; i<rows; i++) {
    for (int j=0; j<cols; j++) {
      fprintf(stdout, "%d  ", mt[i][j]);
      if ((j+1)%8==0)
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }
}

void print_matrix_float(std::vector<std::vector<float>> &mt,
                        int rows,
                        int cols)
{
  for (int i=0; i<rows; i++) {
    for (int j=0; j<cols; j++) {
      fprintf(stdout, "%f  ", mt[i][j]);
      if ((j+1)%8==0)
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }
}

void print_matrix_double(std::vector<std::vector<double>> &mt,
                         int rows,
                         int cols)
{
  for (int i=0; i<rows; i++) {
    for (int j=0; j<cols; j++) {
      fprintf(stdout, "%f  ", mt[i][j]);
      if ((j+1)%8==0)
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
  }
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*---------------    Numerical Programs               ---------------*/
/*---------------    LU decomposition programs        ---------------*/
/*---------------    Calculate matrix inversion       ---------------*/
/*---------------    Calculate matrix determinant     ---------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

// compute 2d matrix determinant
float mat_det_float( std::vector<std::vector<float>> mt,
                    int dim)
{
  int i,j,k,m,n;
  float res;
  
  if (dim==1)
    return(mt[0][0]);
  
  std::vector<std::vector<float>> submt(dim-1, std::vector<float>(dim-1, 0.0));
  
  for (int i=1; i<dim; i++) {
    for (int j=1; j<dim; j++) {
      submt[i-1][j-1] = mt[i][j];
    }
  }
  
  n=1;
  res = 0.0;
  
  for (int i=0; i<dim; i++) {
    res += (n*mt[i][0]*mat_det_float(submt, dim-1));
    n = -n;
    if (i==dim-1)
      continue;
    for (int j=1; j<dim; j++) {
      submt[i][j-1] = mt[i][j];
    }
  }
  
  return(res);
}

double mat_det_double(std::vector<std::vector<double>> &mt,
                      int dim)
{
  double res;
  
  if (dim==1)
    return(mt[0][0]);
  
  int n = dim-1;
  std::vector<std::vector<double>> submt(dim-1, std::vector<double>(dim-1));
  
  for (int i=1; i<dim; i++) {
    for (int j=1; j<dim; j++) {
      submt[i-1][j-1] = mt[i][j];
    }
  }
  
  n=1;
  res = 0.0;
  
  for (int i=0; i<dim; i++) {
    res += (n*mt[i][0]*mat_det_double(submt, dim-1));
    n = -n;
    if (i==dim-1)
      continue;
    for (int j=1; j<dim; j++) {
      submt[i][j-1] = mt[i][j];
    }
  }
  return(res);
}

unsigned char ludcmp_float(std::vector<std::vector<float>> &a,
                           int n,
                           std::vector<int> &indx,
                           float &d)
{
  int imax = 0;
  float big, dum,sum,temp;
  std::vector<float> vv(n, 0.0);
  float TINY=1e-20;
  
  d = 1.0;
  for (int i=0; i<n; i++) {
    big =0.0;
    for (int j=0; j<n; j++)
      if ((temp=fabs(a[i][j]))>big)
        big = temp;
    if (big==0.0) {
      fprintf(stderr, "Singular matrix in ludcmp_float\n");
      return(2); /* 2 stands for singular matrix */
    }
    vv[i]=1.0/big;
  }
  
  for (int j=0; j<n; j++) {
    for (int i=0; i<j; i++) {
      sum = a[i][j];
      for (int k=0; k<i; k++)
        sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (int i=j; i<n; i++) {
      sum = a[i][j];
      for (int k=0; k<j; k++)
        sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum=vv[i]*fabs(sum))>=big) {
        big = dum;
        imax = i;
      }
    }
    if (j!=imax) {
      for (int k=0; k<n; k++) {
        dum = a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j]==0.0) a[j][j] = TINY;
    if (j!=n-1) {
      dum = 1.0/(a[j][j]);
      for (int i=j+1; i<n; i++)
        a[i][j] *= dum;
    }
  }
  return(1);
}


unsigned char ludcmp_double(std::vector<std::vector<double>> &a,
                            int n,
                            std::vector<int> &indx,
                            double &d)
{
  int imax = 0;
  double big, dum,sum,temp;
  std::vector<double> vv(n);
  double TINY=1e-20;
  
  d = 1.0;
  for (int i=0; i<n; i++) {
    big =0.0;
    for (int j=0; j<n; j++)
      if ((temp=fabs(a[i][j]))>big)
        big = temp;
    if (big==0.0) {
      fprintf(stderr, "Singular matrix in ludcmp_double\n");
      return(2); /* 2 stands for singular matrix */
    }
    vv[i]=1.0/big;
  }
  
  for (int j=0; j<n; j++) {
    for (int i=0; i<j; i++) {
      sum = a[i][j];
      for (int k=0; k<i; k++)
        sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (int i=j; i<n; i++) {
      sum = a[i][j];
      for (int k=0; k<j; k++)
        sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum=vv[i]*fabs(sum))>=big) {
        big = dum;
        imax = i;
      }
    }
    if (j!=imax) {
      for (int k=0; k<n; k++) {
        dum = a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j]==0.0) a[j][j] = TINY;
    if (j!=n-1) {
      dum = 1.0/(a[j][j]);
      for (int i=j+1; i<n; i++)
        a[i][j] *= dum;
    }
  }
  return(1);
}

void lubksb_float(std::vector<std::vector<float>> &a,
                  int n,
                  std::vector<int> &indx,
                  std::vector<float> &b)
{
  int ii=-1, ip;
  float sum;
  
  for (int i=0; i<n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (int j=ii; j<i; j++) sum -= a[i][j]*b[j];
    else if (sum != 0.0) ii=i;
    b[i]=sum;
  }
  for (int i=n-1; i>=0; i--) {
    sum=b[i];
    for (int j=i+1; j<n; j++) sum-= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


void lubksb_double(std::vector<std::vector<double>> &a,
                   int n,
                   std::vector<int> &indx,
                   std::vector<double> &b)
{
  int i, ii=-1, ip,j;
  double sum;
  
  for (int i=0; i<n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (int j=ii; j<i; j++) sum -= a[i][j]*b[j];
    else if (sum != 0.0) ii=i;
    b[i]=sum;
  }
  for (int i=n-1; i>=0; i--) {
    sum=b[i];
    for (int j=i+1; j<n; j++) sum-= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


unsigned char mat_inv_float(std::vector<std::vector<float>> &mt,
                            std::vector<std::vector<float>> &y,
                            int dim)
{
  float d;
  std::vector<std::vector<float>> a(dim, std::vector<float>(dim, 0.0));
  
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      a[i][j] = mt[i][j];
    }
  }
  
  std::vector<float> col(dim);
  std::vector<int> indx(dim);
  
  ludcmp_float(a,dim,indx,d);
  for (int j=0; j<dim; j++) {
    for (int i=0; i<dim; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_float(a,dim,indx,col);
    for (int i=0; i<dim;i++)
      y[i][j]=col[i];
  }
  
  return(1);
}

unsigned char mat_inv_double(std::vector<std::vector<float>> &mt,
                             std::vector<std::vector<float>> &y,
                             int dim)
{
  double d;
  
  std::vector<std::vector<double>> a(dim, std::vector<double>(dim,0.0));
  
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      a[i][j] = mt[i][j];
    }
  }
  
  std::vector<double> col(dim);
  std::vector<int> indx(dim);
  
  ludcmp_double(a,dim,indx,d);
  for (int j=0; j<dim; j++) {
    for (int i=0; i<dim; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_double(a,dim,indx,col);
    for (int i=0; i<dim;i++)
      y[i][j]=col[i];
  }
  return(1);
}

float mat_det_ludcmp_float(std::vector<std::vector<float>> &mt,
                           int dim)
{
  float d;
  
  std::vector<std::vector<float>> a(dim, std::vector<float>(dim, 0.0));
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      a[i][j] = mt[i][j];
    }
  }
  std::vector<int> indx(dim, 0);
  
  ludcmp_float(a,dim,indx,d);
  for (int j=0; j<dim; j++) d *= a[j][j];

  return(d);
}

double mat_det_ludcmp_double(std::vector<std::vector<float>> &mt,
                             int dim)
{
  double d;
  
  std::vector<std::vector<double>> a(dim,std::vector<double>(dim, 0.0));
  
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      a[i][j] = mt[i][j];
    }
  }
  std::vector<int> indx(dim);
  
  ludcmp_double(a,dim,indx,d);
  
  for (int j=0; j<dim; j++) d *= a[j][j];

  return(d);
}


/** compute the determinant and at the same provide the inverse matrix **/
/** combining these two operations save computation since they share   **/
/** the call to ludcmp_double                      **/

unsigned char mat_det_inv_double(std::vector<std::vector<double>> &mt,
                                 std::vector<std::vector<double>> &y,
                                 double &det,
                                 int dim)
{
  double d;
  std::vector<std::vector<double>> a(mt);
  
  /** initialize matrix determinant **/
  det=0.0;
  std::vector<double> col(dim, 0.0);
  std::vector<int> indx(dim, 0);
  
  int m=ludcmp_double(a,dim,indx,d);
  if (m==2) {  /** singular matrix **/
    det=0.0;
    return(2);
  }
  
  for (int j=0; j<dim; j++) d *= a[j][j];
  det=d;
  
  for (int j=0; j<dim; j++) {
    for (int i=0; i<dim; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_double(a,dim,indx,col);
    for (int i=0; i<dim;i++)
      y[i][j]=col[i];
  }
  return(1);
}



unsigned char mat_det_inv_float(std::vector<std::vector<float>> &mt,
                                std::vector<std::vector<float>> &y,
                                float &det,
                                int dim)
{
  float d;
  
  /** initialize matrix determinant **/
  det=0.0;
  std::vector<std::vector<float>> a(mt);
  
  std::vector<float> col(dim);
  std::vector<int> indx(dim);
  
  int m=ludcmp_float(a,dim,indx,d);
  if (m==2) {  /** singular matrix **/
    det=0.0;
    return(2);
  }
  
  for (int j=0; j<dim; j++) d *= a[j][j];
  det=d;
  
  for (int j=0; j<dim; j++) {
    for (int i=0; i<dim; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_float(a,dim,indx,col);
    for (int i=0; i<dim;i++)
      y[i][j]=col[i];
  }
  return(1);
}

