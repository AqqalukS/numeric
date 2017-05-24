#include "jacobi.h"
#include <stdlib.h>
#define RND ((double)rand() / RAND_MAX)

int main (int argc, char** argv) {
	int n = atof((argv[1]));
	// Exercise A
	// allocate matrices 
	matrix *A = matrix_alloc (n, n);
	matrix *V = matrix_alloc (n, n);
	vector *b = vector_alloc (n);

	// set values on selected matrces and vectors
	for (int i = 0; i < n; i++) {
		vector_set (b, i, RND);
		for (int j = i; j < n; j++) {
			matrix_set (A, i, j, RND);
			matrix_set (A, j, i, matrix_get (A, i, j));
		}
	}
	int rotations = 0;
	for (int i = 1; i < n; i++) {
		rotations += jacobi_cyclic_by_row (A, b, V, i);
	}
	printf("%i %i \n", n, rotations);
	// free stuff
	matrix_free (V);
	matrix_free (A);
	vector_free (b);
	return 0;
}
