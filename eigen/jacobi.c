#include "jacobi.h"


int jacobi_cyclic (matrix *A, vector *x, matrix *V) {
	assert(A->size1 == A->size2);
	assert(A->size1 == x->size);
	assert(V->size1 == V->size2);
	assert(A->size1 == A->size1);
	
	int changed, rotations = 0, n = A->size1;
	// save 
	for (int i = 0; i < n; i++) {
		vector_set (x, i, matrix_get (A, i, i));
	}

	matrix_set_identity (V);
	do {
		changed = 0; int p, q;
		for (p = 0; p < n; p++) {
			for (q = p+1; q < n; q++) {
				rotations++; 
				double app = vector_get (x, p);
				double aqq = vector_get (x, q);
				double apq = matrix_get (A, p, q);
				double phi = 0.5*atan2(2*apq,aqq-app);
				double c = cos(phi), s = sin(phi);

				double app1 = c*c*app - 2*s*c*apq + s*s*aqq;
				double aqq1 = s*s*app + 2*s*c*apq + c*c*aqq;
				if (app1 != app || aqq1 != aqq) {
					changed = 1;
					vector_set (x, p, app1);
					vector_set (x, q, aqq1);
					matrix_set (A, p, q, 0.0);
					for (int i = 0; i < p; i++) {
						double aip = matrix_get (A, i, p);
						double aiq = matrix_get (A, i, q);
						matrix_set (A, i, p, c*aip - s*aiq);
						matrix_set (A, i, q, c*aiq + s*aip);
					}
					for (int i = p+1; i < q; i++) {
						double api = matrix_get (A, p, i);
						double aiq = matrix_get (A, i, q);
						matrix_set (A, p, i, c*api - s*aiq);
						matrix_set (A, i, q, c*aiq + s*api);
					}
					for (int i = q+1; i < n; i++) {
						double api = matrix_get (A, p, i);
						double aqi = matrix_get (A, q, i);
						matrix_set (A, p, i, c*api - s*aqi);
						matrix_set (A, q, i, c*aqi + s*api);
					}
					for (int i = 0; i < n; i++) {
						double vip = matrix_get (V, i, p);
						double viq = matrix_get (V, i, q);
						matrix_set (V, i, p, c*vip - s*viq);
						matrix_set (V, i, q, c*viq + s*vip);
					}
				}
			}
		}
	} while (changed!=0);
	return rotations;
}

int jacobi_cyclic_by_row (matrix *A, vector *x, matrix *V, int N) {
	assert(A->size1 == A->size2 && 
		A->size1 == x->size &&
		V->size1 == V->size2 &&
		A->size1 == A->size1 &&
		N <= A->size1 &&
		N > 0);
	
	int changed, rotations = 0;
	
	if (N==1) {
		for (int i = 0; i < A->size1; i++) {
			vector_set (x, i, matrix_get (A, i, i));
		}
		matrix_set_identity (V);
	}
	int n = A->size1, p = N-1, q;

	do {
		changed = 0; 
		for (q = p+1; q < n; q++) {
			rotations++;
			double app = vector_get (x, p);
			double aqq = vector_get (x, q);
			double apq = matrix_get (A, p, q);
			double phi = 0.5*atan2(2*apq,aqq-app);
			double c = cos(phi), s = sin(phi);

			double app1 = c*c*app - 2*s*c*apq + s*s*aqq;
			double aqq1 = s*s*app + 2*s*c*apq + c*c*aqq;
			if (app1 != app || aqq1 != aqq) {
				changed = 1;
				vector_set (x, p, app1);
				vector_set (x, q, aqq1);
				matrix_set (A, p, q, 0.0);
				for (int i = 0; i < p; i++) {
					double aip = matrix_get (A, i, p);
					double aiq = matrix_get (A, i, q);
					matrix_set (A, i, p, c*aip - s*aiq);
					matrix_set (A, i, q, c*aiq + s*aip);
				}
				for (int i = p+1; i < q; i++) {
					double api = matrix_get (A, p, i);
					double aiq = matrix_get (A, i, q);
					matrix_set (A, p, i, c*api - s*aiq);
					matrix_set (A, i, q, c*aiq + s*api);
				}
				for (int i = q+1; i < n; i++) {
					double api = matrix_get (A, p, i);
					double aqi = matrix_get (A, q, i);
					matrix_set (A, p, i, c*api - s*aqi);
					matrix_set (A, q, i, c*aqi + s*api);
				}
				for (int i = 0; i < n; i++) {
					double vip = matrix_get (V, i, p);
					double viq = matrix_get (V, i, q);
					matrix_set (V, i, p, c*vip - s*viq);
					matrix_set (V, i, q, c*viq + s*vip);
				}
			}
		}
	} while (changed!=0);
	return rotations;
}

int jacobi_classic (matrix *A, vector *x, matrix *V) {
	assert(A->size1 == A->size2 && 
		A->size1 == x->size &&
		V->size1 == V->size2 &&
		A->size1 == A->size1);
	
	int changed, rotations = 0, n = A->size1;
	// save 
	for (int i = 0; i < n; i++) {
		vector_set (x, i, matrix_get (A, i, i));
	}

	matrix_set_identity (V);
	do {
		changed = 0; int i, j, p = 0, q = 0;
		for (i = 0; i < n; i++) {
			for (j = i+1; j < n; j++) {
				if (matrix_get (A, i, j) > matrix_get (A, p, q)) {
					p = i; q = j;
				}
			}
		}
		rotations++; 
		double app = vector_get (x, p);
		double aqq = vector_get (x, q);
		double apq = matrix_get (A, p, q);
		double phi = 0.5*atan2(2*apq,aqq-app);
		double c = cos(phi), s = sin(phi);

		double app1 = c*c*app - 2*s*c*apq + s*s*aqq;
		double aqq1 = s*s*app + 2*s*c*apq + c*c*aqq;
		if (app1 != app || aqq1 != aqq) {
			changed = 1;
			vector_set (x, p, app1);
			vector_set (x, q, aqq1);
			matrix_set (A, p, q, 0.0);
			for (int i = 0; i < p; i++) {
				double aip = matrix_get (A, i, p);
				double aiq = matrix_get (A, i, q);
				matrix_set (A, i, p, c*aip - s*aiq);
				matrix_set (A, i, q, c*aiq + s*aip);
			}
			for (int i = p+1; i < q; i++) {
				double api = matrix_get (A, p, i);
				double aiq = matrix_get (A, i, q);
				matrix_set (A, p, i, c*api - s*aiq);
				matrix_set (A, i, q, c*aiq + s*api);
			}
			for (int i = q+1; i < n; i++) {
				double api = matrix_get (A, p, i);
				double aqi = matrix_get (A, q, i);
				matrix_set (A, p, i, c*api - s*aqi);
				matrix_set (A, q, i, c*aqi + s*api);
			}
			for (int i = 0; i < n; i++) {
				double vip = matrix_get (V, i, p);
				double viq = matrix_get (V, i, q);
				matrix_set (V, i, p, c*vip - s*viq);
				matrix_set (V, i, q, c*viq + s*vip);
			}
		}
	} while (changed!=0);
	return rotations;
}
