#ifndef HAVE_LINEQ_H

#include <math.h>
#include <matrix.h>
#include <vector.h>

void backsub (matrix *U, vector *v);
void qr_gs_decomp (matrix *A, matrix *R);
void qr_gs_solve (const matrix *Q, const matrix *R, vector *b, vector *x);
void qr_gs_inverse (const matrix *Q, const matrix *R, matrix *B);

#define HAVE_LINEQ_H 
#endif
