#ifndef INTERP_H
#include "vector.h"

typedef struct {int n; double *x, *y, *b, *c;} qspline;
typedef struct {int n; double *x, *y, *b, *c, *d;} cspline;

// Linear 
double linterp (vector *x, vector *y, double z);

// Quadratic spline
qspline * qspline_alloc (vector *x, vector *y);
double qspline_evaluate (qspline *s, double z);
void qspline_free (qspline *s);

// Cubic spline
cspline * cspline_alloc (vector *x, vector *y);
double cspline_evaluate (cspline *s, double z);
void cspline_free (cspline *s);
#define INTERP_H
#endif
