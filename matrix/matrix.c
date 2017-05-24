#include "matrix.h"

matrix * matrix_alloc (int n, int m) {
	matrix *A = (matrix*) malloc(sizeof(matrix));
	A->size1 = n; A->size2 = m;
	A->data = (double*) malloc(n*m*sizeof(double));
	return A;
}


void matrix_set (matrix *A, int i, int j, double x) {
	A->data[i+j*A->size1] = x;
}

double matrix_get (matrix *A, int i, int j) {
	return A->data[i+j*A->size1];
}

void matrix_free (matrix *A) {
	free(A->data); free(A);
}

double matrix_cv_innerp (matrix *A, int i, matrix *B, int j) {
	assert(A->size1 == A->size1);
	// column vector inner product
	double x = 0;
	for (int k = 0; k < A->size1; k++) {
		x += matrix_get (A, k, i) * matrix_get (B, k, j);
	}
	return x;
}

void matrix_print (matrix *A) {
	int i, j;
	for (i = 0; i < A->size1; i++) {
		for (j = 0; j < A->size2; j++) {
			printf("%7.4f\t", matrix_get (A, i, j));
		}
		printf("\n");
	}
	printf("\n");
}

void matrix_x_matrix (matrix *A, matrix *B, matrix *C) {
	assert (A->size2 == B->size1 && 
		C->size1 == A-> size1 && 
		C->size2 == B->size2);
	for (int i = 0; i < A->size1; i++) {
		for (int j = 0; j < B->size2; j++) {
			double s = 0;
			for (int k = 0; k < A->size2; k++) {
				s += matrix_get (A, i, k) 
						* matrix_get (B, k, j);
			}
			matrix_set (C, i, j, s);
		}
	}
}

void matrixT_x_matrix (matrix *A, matrix *B, matrix *C) {
	assert(A->size1 == B->size1 &&
		C->size1 == A->size2 &&
		C->size2 == B->size2);
	for (int i = 0; i < A->size2; i++) {
		for (int j = 0; j < B->size1; j++) {
			double s = 0;
			for (int k = 0; k < A->size1; k++) {
				s += matrix_get (A, k, i) * matrix_get (B, k, j);
			}
			matrix_set (C, i, j, s);
		}
	}
}


void matrix_x_matrixT (matrix *A, matrix *B, matrix *C) {
	assert(A->size2 == B->size2 &&
		C->size1 == A->size1 &&
		C->size2 == B->size1);
	for (int i = 0; i < C->size1; i++) {
		for (int j = 0; j < C->size2; j++) {
			double s = 0;
			for (int k = 0; k < A->size2; k++) {
				s+=matrix_get (A, i, k) * matrix_get (B, j, k);
			}
			matrix_set (C, i, j, s);
		}
	}
}

void matrix_x_vector (matrix *A, vector *v, vector *b) {
	assert(A->size2 == v->size && A->size1 == b->size);
	for (int i = 0; i < A->size1; i++) {
		double s = 0;
		for (int j = 0; j < A->size2; j++) {
			s += matrix_get (A, i, j) * vector_get (v, j);
		}
		vector_set (b, i, s);
	}
}

void matrixT_x_vector (matrix *A, vector *v, vector *b) {
	assert(A->size1 == v->size && A->size2 == b->size);
	for (int i = 0; i < A->size2; i++) {
		double s = 0;
		for (int j = 0; j < A->size1; j++) {
			s += matrix_get (A, j, i) * vector_get (v, j);
		}
		vector_set (b, i, s);
	}
}

void matrix_set_column (matrix *A, int i, vector *v) {
	assert(A->size1 == v->size);
	for (int j = 0; j < A->size1; j++) {
		matrix_set (A, j, i, vector_get (v, j));
	}
}

void matrix_get_column (matrix *A, int i, vector *v) {
	assert(A->size1 == v->size);
	for (int j = 0; j < A->size1; j++) {
		vector_set (v, j, matrix_get (A, j, i));
	}
}

void matrix_set_identity (matrix *I) {
	assert(I->size1 == I->size2);
	for (int i = 0; i < I->size1; i++) {
		matrix_set (I, i, i, 1.0);
		for (int j = i + 1; j < I->size1; j++) {
			matrix_set (I, i, j, 0);
			matrix_set (I, j, i, 0);
		}
	}
}

