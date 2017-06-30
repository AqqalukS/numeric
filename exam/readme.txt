Lanczos tridiagonal algorithm (Arnoldi iteration)

Exercise A
Implementation of Lanczos algorithm for real symmetric matrices is found in lanczos.c, where the input is a real symmetric matrix A size nxn, vector u of size n. output is tridiagonal real symmetric matrix T size of mxm, and nxm matrix V with orthonormal columns.

	lanczos (matrix *A, vector *u, matrix *V, matrix *T);

The vector v is normalized in the function.

In mainA.c the Lanczos algorithm is used, for simplicity is m = n, so that A = V*T*V^T can be checked.

In the outA.txt, it can be seen that the matrix T is tridiagonal, V^T*V = V*V^T = I indentity and that A = V*T*V^T is true. 

Exercise B
Implementation of Arnoldi iteration for real square matrices is found in arnoldi.c, where the input is a real square matrix A with a size of nxn, vector q size of n. output is an upper Hessenberg matrix H nxn, and a matrix Q nxn which holds sequence of orthonormal column vectors called Arnoldi vectors.
	arnoldi (matrix *A, vector *q, matrix *Q, matrix *H);

The vector q is q_1, and is normalized in the function. 

In mainB.c the Arnoldi iteration is used, and notice that if the input matrix is a real symmetric matrix A, then the output matrix H is the same as the output matrix T in Lanczos algorithm (see outB.txt). Q*H*Q^T = A i checked in outB.txt.
