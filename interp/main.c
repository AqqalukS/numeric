#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interp.h"
#include "vector.h"
int main(int argc, const char *argv[]) {
	int n = atof(argv[1]);

	vector *x = vector_alloc (n);
	vector *y = vector_alloc (n);

	fprintf(stderr, "# datapoints\n");
	for (int i = 0; i < n; i++) {
		vector_set (x, i, M_PI/n*i);
		vector_set (y, i, cos(vector_get(x, i)));
		fprintf(stderr, "%g %g\n", 
				vector_get(x, i), 
				vector_get(y, i));
	}
	
	qspline *qs = qspline_alloc (x, y);
	cspline *cs = cspline_alloc (x, y);
	double dz = 0.001;
	printf("# interpolation \n");
	printf("# x linterp qspline cspline\n");
	for (double z = vector_get(x, 0); z < vector_get(x, n-1); z += dz) {
		double l1 = linterp (x, y, z);
		double l2 = qspline_evaluate (qs, z);
		double l3 = cspline_evaluate (cs, z);
		printf("%g %g %g %g\n", z, l1, l2, l3);
	}
	qspline_free (qs);
	cspline_free (cs);
	vector_free(x);
	vector_free(y);
	return 0;
}
