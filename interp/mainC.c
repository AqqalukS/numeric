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
	
	cspline *cs = cspline_alloc (x, y);
	double dz = 0.001;
	printf("# interpolation \n");
	printf("# x qspline integ deriv\n");
	for (double z = vector_get(x, 0); z < vector_get(x, n-1); z += dz) {
		double l3 = cspline_evaluate (cs, z); 		// Quatratic spline
		double i3 = cspline_integ (cs, z);
		double d3 = cspline_deriv (cs, z);
		printf("%g %g %g %g\n", z, l3, i3, d3);
	}
	cspline_free (cs);
	vector_free(x);
	vector_free(y);
	return 0;
}
