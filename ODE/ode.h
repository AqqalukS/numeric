#ifndef HAVE_ODE_H
#include "vector.h"
void rkstep23 (double t, double h, vector *y, void f(double t, vector *y, vector *dydt), vector *yh, vector *err);
#define HAVE_ODE_H 
#endif
