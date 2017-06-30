#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interp.h"
#include "vector.h"
#define RND ((double)rand() / RAND_MAX)

int main(int argc, const char *argv[]) {
	int n = atof(argv[1]);

	vector *x = vector_alloc (n);
	vector *y = vector_alloc (n);

	fprintf(stderr, "# datapoints\n");
	for (int i = 0; i < n; i++) {
		vector_set (x, i, i * M_PI/20);
		vector_set (y, i, cos(vector_get(x, i)));
		fprintf(stderr, "%g %g\n", 
				vector_get(x, i), 
				vector_get(y, i));
	}
	
	qspline *qs = qspline_alloc (x, y);
	double dz = 0.001;
	printf("# interpolation \n");
	printf("# x qspline integ deriv\n");
	for (double z = vector_get(x, 0); z < vector_get(x, n-1); z += dz) {
		double l2 = qspline_evaluate (qs, z); 		// Quatratic spline
		double i2 = qspline_integ (qs, z);
		double d2 = qspline_deriv (qs, z);
		printf("%g %g %g %g\n", z, l2, i2, d2);
	}
	qspline_free (qs);
	vector_free(x);
	vector_free(y);
	return 0;
}
