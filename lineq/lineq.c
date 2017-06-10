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
	assert(Q->size1 == b->size);
	assert(Q->size2 == R->size1);
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

void givens_decomp (matrix *A) {
	for (int q = 0; q < A->size2; q++) {
		for (int p = q + 1; p < A->size1; p++) {
			double theta = atan2 (matrix_get (A, p, q),
					matrix_get (A, q, q));
			for (int k = q; k < A->size2; k++) {
				double xq = matrix_get (A, q, k);
				double xp = matrix_get (A, p, k);
				matrix_set (A, q, k, xq * cos(theta) + xp * sin(theta));
				matrix_set (A, p, k,-xq * sin(theta) + xp * cos(theta));
			}
			matrix_set (A, p, q, theta);
		}
	}
}

void givens_QTvec (matrix *QR, vector *v) {
	for (int q = 0; q < QR->size2; q++) {
		for (int p = q + 1; p < QR->size1; p++) {
			double theta = matrix_get (QR, p, q);
			double vq = vector_get (v, q);
			double vp = vector_get (v, p);
			vector_set (v, q, vq * cos(theta) + vp * sin(theta));
			vector_set (v, p,-vq * sin(theta) + vp * cos(theta));
		}
	}
}

void givens_solve (matrix *QR, vector *b) {
	givens_QTvec (QR, b);
	backsub (QR, b);
}

void givens_unpack_Q (matrix *QR, matrix *Q) {
	vector *ei = vector_alloc (QR->size1);
	for (int i = 0; i < QR->size1; i++) {
		for (int j = 0; j < ei->size; j++) {
			if (j == i) vector_set (ei, j, 1.);
			else vector_set (ei, j, 0.);
		}
		givens_QTvec (QR, ei);
		for (int j = 0; j < QR->size2; j++) {
			matrix_set (Q, i, j, vector_get (ei, j));
		}
	}
	vector_free (ei);
}

void givens_inverse (matrix *QR, matrix *B) {
	assert(QR->size1 == QR->size2 && QR->size1 == B->size1 
			&& B->size1 == B->size2);

	matrix_set_identity (B);
	vector *b = vector_alloc (QR->size1);
	for (int i = 0; i < QR->size2; i++) {
		matrix_get_column (B, i, b);
		givens_solve (QR, b);
		matrix_set_column (B, i, b);
	}
	vector_free (b);
}
