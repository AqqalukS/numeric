#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ode.h"

void rkstep23 (double t, double h, vector *y, void f(double t, vector *y, vector *dydt), vector *yh, vector *err);
void driver (
	double *t, double b, double *h, vector *y, double acc, double eps, 
	void stepper (
		double t, double h, vector *y, 
		void f(double t, vector *y, vector *dydt), vector *yh, vector *err
		), 
	void f (double t, double *y, double *dydt));

int main (int argc, const char *argv[]) {
	double b = M_PI;
	return 0;
}
