#ifndef __hmm_aw_matrix_h
#define __hmm_aw_matrix_h
#include <vector>
/*-------------------------- mat_simple.c --------------------------*/

/*-------------------------------------------------------------------*/
/*------------------------ Print out Vector -------------------------*/
/*-------------------------------------------------------------------*/
extern void print_vector_uchar(std::vector<unsigned char>& vt);
extern void print_vector_int(std::vector<int>& vt);
extern void print_vector_float(std::vector<float>& vt);
extern void print_vector_double(std::vector<double>& vt);
/*-------------------------------------------------------------------*/
/*------------------------ Print out Matrix -------------------------*/
/*-------------------------------------------------------------------*/
extern void print_matrix_uchar(std::vector<std::vector<unsigned char>> &mt);
extern void print_matrix_int(std::vector<std::vector<int>> &mt);
extern void print_matrix_float(std::vector<std::vector<float>> &mt);
extern void print_matrix_double(std::vector<std::vector<double>> &mt);


/*-------------------------------------------------------------------*/
/*---------------      Numerical Programs             ---------------*/
/*---------------      LU decomposition programs      ---------------*/
/*---------------      Calculate matrix inversion     ---------------*/
/*---------------      Calculate matrix determinant   ---------------*/
/*-------------------------------------------------------------------*/
extern float mat_det_float(std::vector<std::vector<float>> mt,
                           int dim);
extern double mat_det_double(std::vector<std::vector<double>> &mt,
                             int dim);
extern unsigned char ludcmp_float(std::vector<std::vector<float>> &a,
                                  int n,
                                  std::vector<int> &indx,
                                  float &d);
extern unsigned char ludcmp_double(std::vector<std::vector<double>> &a,
                                   int n,
                                   std::vector<int> &indx,
                                   double &d);
extern void lubksb_float(std::vector<std::vector<float>> &a,
                         int n,
                         std::vector<int> &indx,
                         std::vector<float> &b);
extern void lubksb_double(std::vector<std::vector<double>> &a,
                          int n,
                          std::vector<int> &indx,
                          std::vector<double> &b);
extern unsigned char mat_inv_float(std::vector<std::vector<float>> &mt,
                                   std::vector<std::vector<float>> &y,
                                   int dim);
extern unsigned char mat_inv_double(std::vector<std::vector<float>> &mt,
                                    std::vector<std::vector<float>> &y,
                                    int dim);
extern float mat_det_ludcmp_float(std::vector<std::vector<float>> &mt,
                                  int dim);
extern double mat_det_ludcmp_double(std::vector<std::vector<float>> &mt,
                                    int dim);
unsigned char mat_det_inv_double(std::vector<std::vector<double>> &mt,
                                 std::vector<std::vector<double>> &y,
                                 double &det,
                                 int dim);
unsigned char mat_det_inv_float(std::vector<std::vector<float>> &mt,
                                std::vector<std::vector<float>> &y,
                                float &det,
                                int dim);
#endif
