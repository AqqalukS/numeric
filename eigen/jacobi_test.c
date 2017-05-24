#include "jacobi.h"
#include <stdlib.h>
#define RND ((double)rand() / RAND_MAX)

int main (int argc, char** argv) {
	int n = atof(argv[1]);
	// Exercise A
	// allocate matrices 
	matrix *A = matrix_alloc (n, n);
	matrix *V = matrix_alloc (n, n);
	matrix *R1 = matrix_alloc (n, n);
	matrix *R2 = matrix_alloc (n, n);
	matrix *R3 = matrix_alloc (n, n);
	matrix *R4 = matrix_alloc (n, n);
	vector *b = vector_alloc (n);

	// set values on selected matrces and vectors
	for (int i = 0; i < n; i++) {
		vector_set (b, i, RND);
		for (int j = i; j < n; j++) {
			matrix_set (A, i, j, RND);
			matrix_set (A, j, i, matrix_get (A, i, j));
		}
	}
	printf("# Exercise A\n");

	printf("# Symetric square matrix A\n");
	matrix_print (A);

	// Jacobi diagonalization
	int *jacobi = jacobi_cyclic (A, b, V);
	
	printf("# Distroyed upper triangular matrix A\n");
	matrix_print (A);

	// restore upper triangular of A
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			matrix_set (A, i, j, matrix_get (A, j, i));
		}
	}
	printf("# Number of sweeps = %i\n", jacobi[0]);
	printf("# Number of rotations = %i\n", jacobi[1]);

	printf("# Orthogonal matrix V of eigenvectors \n");
	matrix_print (V);

	matrixT_x_matrix (V, V, R1);
	printf("# V^T times V \n");
	matrix_print (R1);

	matrix_x_matrixT (V, V, R2);
	printf("# V times V^T \n");
	matrix_print (R2);

	matrixT_x_matrix (V, A, R3);
	matrix_x_matrix (R3, V, R4);
	printf("# V^T*A*V = D \n");
	matrix_print (R4);

	printf("# Eigenvalues \n");
	vector_print (b);

	// Exercise B
	int rotations;
	printf("# Exercise B\n");
	printf("# Sweeps row by row\n");
	for (int i = 1; i < 3; i++) {
		rotations = jacobi_cyclic_by_row (A, b, V, i);
		printf("rotations row(%i) = %i\n", i, rotations);
		 matrix_print (A);
	}

	printf("# The first row \n");
	printf("# To calculate the largest eigenvalue first, you should add pi/2 to phi\n");
	printf("# matrix V\n");
	matrix_print (V);

	matrixT_x_matrix (V, V, R1);
	printf("# V^T times V \n");
	matrix_print (R1);

	matrix_x_matrixT (V, V, R2);
	printf("# V times V^T \n");
	matrix_print (R2);

	// restore upper triangular of A
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			matrix_set (A, i, j, matrix_get (A, j, i));
		}
	}

	matrixT_x_matrix (V, A, R3);
	matrix_x_matrix (R3, V, R4);
	printf("# V^T*A*V = D \n");
	matrix_print (R4);

	vector_print (b);

	fprintf(stderr,"\n");
	// free stuff
	matrix_free (V);
	matrix_free (A);
	matrix_free (R1);
	matrix_free (R2);
	matrix_free (R3);
	matrix_free (R4);
	vector_free (b);
	free (jacobi);
	return 0;
}
