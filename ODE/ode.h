#ifndef HAVE_ODE_H
#include "vector.h"
void rkstep23 (double t, double h, vector *y, void f(double t, vector *y, vector *dydt), vector *yh, vector *err);
void driver (
	double t, double b, double h, vector *y, double acc, double eps, 
	void stepper (
		double t, double h, vector *y, 
		void f(double t, vector *y, vector *dydt), vector *yh, vector *err
		), 
	void f (double t, vector *y, vector *dydt));
#define HAVE_ODE_H 
#endif
