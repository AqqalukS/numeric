#include <math.h>
#include "vector.h"
#include "matrix.h"

void driver (
	double *t, double b, double *h, vector *y, double acc, double eps, 
	void stepper (
		double t, double h, vector *y, 
		void f(double t, vector *y, vector *dydt), vector *yh, vector *err
		), 
	void f (double t, double *y, double *dydt)) {
	int n = y->size, k = 0;
	double tau, e, s;
	vector *yh = vector_alloc (n);
	vector *err = vector_alloc (n);
	while (t < b) {
		:
