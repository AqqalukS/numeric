#include "lineq.h"

void backsub (matrix *U, vector *v) {
	assert(U->size1 == v->size);
	for (int i = v->size - 1; i >= 0; i--) {
		double s = 0;
		for (int k = i + 1; k < v->size; k++) {
			s += matrix_get (U, i, k) * vector_get (v, k);
		}
		vector_set (v, i, (vector_get(v, i)-s)/matrix_get (U, i, i));
	}
}
void qr_gs_decomp (matrix *A, matrix *R) {
	assert(A->size2 == R->size1 && R->size1 == R->size2);
	int i, j, k, n = A->size1, m = A->size2;
	for (i = 0; i < m; i++) { 
		double aa = 0;
		for (k = 0; k < n; k++) {
			aa += matrix_get (A, k, i) * matrix_get (A, k, i);
		}
		matrix_set (R, i, i, sqrt(aa));
		for (k = 0; k < n; k++) {
			matrix_set (A, k, i, matrix_get (A, k, i) 
					/ matrix_get (R, i, i));
		}	
		for (j = i + 1; j < m; j++) {
			double qa = 0; 
			for (k = 0; k < n; k++) {
				qa += matrix_get (A, k, i) * matrix_get (A, k, j);
			}
			matrix_set (R, i, j, qa);
			for (k = 0; k < n; k++) {
				matrix_set (A, k, j, matrix_get (A, k, j) - 
						matrix_get (A, k, i) * qa);
			}
		}
	}
}

void qr_gs_solve (const matrix *Q, const matrix *R, vector *b, vector *x) {
	assert(Q->size1 == b->size &&
		Q->size2 == R->size1);
	matrixT_x_vector (Q, b, x);
	backsub (R, x);
}

void qr_gs_inverse (const matrix *Q, const matrix *R, matrix *B) {
	assert(Q->size1 == Q->size2 && Q->size1 == B->size1 
			&& B->size1 == B->size2);

	matrix_set_identity (B);

	vector *x = vector_alloc (Q->size1);
	vector *b = vector_alloc (Q->size1);
	for (int i = 0; i < Q->size2; i++) {
		matrix_get_column (B, i, b);
		qr_gs_solve (Q, R, b, x);
		matrix_set_column (B, i, x);
	}
	vector_free (x);
	vector_free (b);
}
