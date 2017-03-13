#ifndef __hmm_aw_matrix_h
#define __hmm_aw_matrix_h
/*-------------------------- mat_simple.c --------------------------*/

/*------------------------------------------------------------------*/
/*-------------------------- General vector ------------------------*/
/*------------------------------------------------------------------*/

extern unsigned char * vector_uchar( int );
extern float * vector_float( int );
extern double * vector_double( int );
extern int * vector_int( int);

/*------------------------------------------------------------------*/
/*---------------- General 2 dimension matrix ----------------------*/
/*------------------------------------------------------------------*/

extern unsigned char ** matrix_2d_uchar( int , int );
extern float ** matrix_2d_float( int , int );
extern double ** matrix_2d_double( int , int );
extern int ** matrix_2d_int( int , int );

/*--------------------------------------------------------------------*/
/*------------ Generate 3 dimension matrix ---------------------------*/
/*--------------------------------------------------------------------*/

extern unsigned char *** matrix_3d_uchar( int ,int , int );
extern float *** matrix_3d_float( int ,int , int );
extern double *** matrix_3d_double( int ,int , int );
extern int *** matrix_3d_int( int , int , int );

/*-------------------------------------------------------------------*/
/*--------------- Free 2 dimension matrix ---------------------------*/
/*-------------------------------------------------------------------*/
extern void free_matrix_2d_uchar(unsigned char **, int );
extern void free_matrix_2d_float(float **, int );
extern void free_matrix_2d_double(double **, int );
extern void free_matrix_2d_int(int **, int );

/*-------------------------------------------------------------------*/
/*--------------- Free 3 dimension matrix ---------------------------*/
/*-------------------------------------------------------------------*/
extern void free_matrix_3d_uchar(unsigned char ***, int , int );
extern void free_matrix_3d_float(float ***, int , int );
extern void free_matrix_3d_double(double ***, int , int );
extern void free_matrix_3d_int(int ***, int , int );

/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for Vector ---------------------*/
/*-------------------------------------------------------------------*/
extern void memcpy_1d_uchar(unsigned char *, int , unsigned char );
extern void memcpy_1d_int(int *, int , int );
extern void memcpy_1d_float(float *, int , float );
extern void memcpy_1d_double(double *, int , double );

/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for 2D matrix ------------------*/
/*-------------------------------------------------------------------*/
extern void memcpy_2d_uchar(unsigned char **, int , int , unsigned char );
extern void memcpy_2d_int(int **, int , int , int );
extern void memcpy_2d_float(float **, int , int , float );
extern void memcpy_2d_double(double **, int , int , double );

/*-------------------------------------------------------------------*/
/*----------------------- Set Memory for 3D matrix ------------------*/
/*-------------------------------------------------------------------*/
extern void memcpy_3d_uchar(unsigned char ***, int, int, int,unsigned char);
extern void memcpy_3d_int(int ***, int , int , int , int);
extern void memcpy_3d_float(float ***, int , int , int , float);
extern void memcpy_3d_double(double ***, int , int , int , double);

/*-------------------------------------------------------------------*/
/*------------------------ Vector copy ------------------------------*/
/*-------------------------------------------------------------------*/
extern void vector_cpy_uchar(unsigned char *, unsigned char *, int );
extern void vector_cpy_int(int *, int *, int );
extern void vector_cpy_float(float *, float *, int );
extern void vector_cpy_double(double *, double *, int );

/*-------------------------------------------------------------------*/
/*------------------------ Matrix copy ------------------------------*/
/*-------------------------------------------------------------------*/
extern void matrix_2d_cpy_uchar(unsigned char **, unsigned char **, 
		      int , int );
extern void matrix_2d_cpy_int(int **, int **, int , int );
extern void matrix_2d_cpy_float(float **, float **, int , int );
extern void matrix_2d_cpy_double(double **, double **, int , int );

/*-------------------------------------------------------------------*/
/*------------------------ Print out Matrix -------------------------*/
/*-------------------------------------------------------------------*/
extern void print_matrix_uchar(unsigned char **, int , int );
extern void print_matrix_int(int **, int , int );
extern void print_matrix_float(float **, int , int );
extern void print_matrix_double(double **, int , int );

/*-------------------------------------------------------------------*/
/*---------------      Numerical Programs             ---------------*/
/*---------------    LU decomposition programs        ---------------*/
/*---------------      Calculate matrix inversion     ---------------*/
/*---------------      Calculate matrix determinant   ---------------*/
/*-------------------------------------------------------------------*/
extern float mat_det_float(float **mt, int dim);
extern double mat_det_double(double **mt, int dim);
extern unsigned char ludcmp_float(float **, int , int *, float *);
extern unsigned char ludcmp_double(double **, int , int *, double *);
extern void lubksb_float(float **, int , int *, float *);
extern void lubksb_double(double **, int , int *, double *);
extern unsigned char mat_inv_float(float **, float **, int );
extern unsigned char mat_inv_double(double **, double **, int );
extern float mat_det_ludcmp_float(float **, int );
extern double mat_det_ludcmp_double(double **, int );
unsigned char mat_det_inv_double(double **mt, double **y, double *det,
                                 int dim);
#endif
