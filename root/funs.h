#ifndef HAVE_FUNS_H
#include <math.h>
#include "vector.h"
#include "matrix.h"
void f1 (vector *v, vector *df);
void J_f1 (vector *v, matrix *J);
void rosen (vector *v, vector *df);
void himmel (vector *v, vector *df);
#define HAVE_FUNS_H 
#endif
