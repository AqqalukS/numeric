#ifndef HAVE_FUNS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
double funs (int i, double x);
double fit (double x, double f(int, double), vector *c);
double df (double x, double f(int, double), matrix *S);
#define HAVE_FUNS_H 
#endif
