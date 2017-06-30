#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "ode.h"
#define stepper rkstep23

int main (int argc, const char *argv[]) {
	void f (double t, vector *y, vector *dydt) {
		vector_set (dydt, 0,  vector_get (y, 1));
		vector_set (dydt, 1, -vector_get (y, 0));
	}

	double b = 8*M_PI, t = 0, h = 0.1, acc = 1e-6, eps = 1e-6, dt = M_PI/4;
	int n = 2;
	vector *y = vector_alloc (n);
	vector_set (y, 0, 0);
	vector_set (y, 1, 1);
	fprintf(stderr, "h = %g, t = %g\n", h, t);
	driver (t, b, h, y, acc, eps, stepper, f);
	vector_print (y);

	return 0;
}
