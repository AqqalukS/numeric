#ifndef HAVE_JACOBI_H
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"

int *jacobi_cyclic (matrix *A, vector *e, matrix *V);
int jacobi_cyclic_by_row (matrix *A, vector *e, matrix *V, int N);
int jacobi_classic (matrix *A, vector *e, matrix *V);
#define HAVE_JACOBI_H
#endif
