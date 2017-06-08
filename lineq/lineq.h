#ifndef HAVE_LINEQ_H

#include <math.h>
#include <matrix.h>
#include <vector.h>

void backsub (matrix *U, vector *v);
void qr_gs_decomp (matrix *A, matrix *R);
void qr_gs_solve (const matrix *Q, const matrix *R, vector *b, vector *x);
void qr_gs_inverse (const matrix *Q, const matrix *R, matrix *B);
void givens_decomp (matrix *A);
void givens_QTvec (matrix *QR, vector *v);
void givens_solve (matrix *QR, vector *b);
void givens_unpack_Q (matrix *QR, matrix *Q);

#define HAVE_LINEQ_H 
#endif
