#ifndef HAVE_ROOT_H

#include "vector.h"
#include "matrix.h"

int newton_root (void f(vector *, vector *), vector *x, double dx, double eps);
int newton_jacobian (void f(vector *, vector *), void jacobian(vector *, matrix *), vector *x,				double epsilon);

#define HAVE_ROOT_H 
#endif
