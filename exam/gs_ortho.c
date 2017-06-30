#include <math.h>
#include <assert.h>
#include "vector.h"
#include "matrix.h"

void orthonormalize (vector *v, matrix *V, int m) {
	int n = V->size1;
	assert (m > n);
	assert (v->size == V->size1);

	vector *ui = vector_alloc (n);
	vector *uj = vector_alloc (n);
	matrix *U = matrix_alloc (n, m);
	double UjUj, UiUj;

	matrix_set_column (V, m, v);

	matrix_get_column (V, 0, ui);
	
	UjUj = vectorT_x_vector (ui, ui);
	vector_scale (ui, 1/sqrt(UjUj));

	matrix_set_column (U, 0, ui);

	for (int i = 1; i < m; i++) {
		matrix_get_column (V, i, ui);
		for (int j = 0; j < i-1; j++) {
			matrix_get_column (U, j, uj);
			UiUj = vectorT_x_vector (ui, uj);
			UjUj = vectorT_x_vector (uj, uj);
			vector_scale (uj, UiUj/UjUj);
			vector_minus_vector (ui, uj);
		}
		UjUj = vectorT_x_vector (ui, ui);
		vector_scale (ui, 1/sqrt(UjUj));
		matrix_set_column (U, i, ui);
	}
	matrix_get_column (U, m, v);

	matrix_free (U);
	vector_free (ui);
	vector_free (uj);
}
