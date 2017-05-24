#ifndef HAVE_MATRIX_H
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector.h>

typedef struct {int size1, size2; double *data;} matrix;

matrix * matrix_alloc (int n, int m);
void matrix_set (matrix *A, int i, int j, double x);
double matrix_get (matrix *A, int i, int j);
void matrix_free (matrix *A);
double matrix_cv_innerp (matrix *A, int i, matrix *B, int j);
void matrix_print (matrix *A);
void matrix_x_matrix (matrix *A, matrix *B, matrix *C);
void matrixT_x_matrix (matrix *A, matrix *B, matrix *C);
void matrix_x_matrixT (matrix *A, matrix *B, matrix *C);
void matrix_x_vector (matrix *A, vector *v, vector *b);
void matrixT_x_vector (matrix *A, vector *v, vector *b);
void matrix_set_column (matrix *A, int i, vector *v);
void matrix_get_column (matrix *A, int i, vector *v);
void matrix_set_identity (matrix *I);

#define HAVE_MATRIX_H 
#endif
