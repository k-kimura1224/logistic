#ifndef MYVECTOR_H
#define MYVECTOR_H


/******************************************************
Functions:
	Com_LS_dposv: compute a linear system
	Com_LS_dgesv: compute a linear system
	Com_mat_AtA: compute A^T A
	Com_scal: compute y := (alpha)x
	Com_linecomb: compute z := (alpha)x + (beta)y
	Com_mat_Ax: compute z := A x
	Com_mat_Atx: compute z := A^T x
	Com_det: compute determinant matrix
	Com_eigen: compute eigenvalues and eigenvectors
	Equal: return true if give values are equal
Template Functions:
	TraMat: transpose matrix
	printv: output a vector 
	printM: output a matrix
	Gen_ZeroVec: generate a zero vector
	count: count positive values
	Bubble_Sort: run the bubble sort
Inline Functions:
	Com_nrm: compute 2-norm
	Copy_vec: copy a vector
	Com_dot: compute a inner product
******************************************************/

#if defined(__APPLE__)
#include <cblas.h>
#else
extern "C" {
	void cblas_dcopy( int n, const double *x, int incX, double *y,  int incY); 
}
extern "C" {
	double cblas_dnrm2( int n, const double *x, int incx);
}
extern "C" {
	double cblas_ddot( int n, const double *x, int incx, const double *y, int incy);
}
#endif
// Functions:

/* compute Ax=b with Cholesky decomposition */
/* A: positive definite symmetric matrix */
/* if successful, return 0 */
int Com_LS_dposv(
	const double	*A,		/* A[n*n] */
	const double	*b,		/* b[n]	*/
	int		n,					/* size */
	double			*x			/* solution */
	);

// compute Ax=b with LU decomposition
int Com_LS_dgesv(
	const double	*A,		/* A[n*n] */
	const double	*b,		/* b[n]	*/
	int		n,					/* size */
	double			*x			/* solution */
	);

// compute Ax=b
// A: symmetric matrix
int Com_LS_dsysv(
	const double	*A,		/* A[n*n] */
	const double	*b,		/* b[n]	*/
	int		n,					/* size */
	double			*x			/* solution */
	);

/* y := (alpha)x */
void Com_scal(
	const double	*x,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	double			*y
	);

// compute linear combination
// z := (alpha)x + (beta)y 
void Com_linecomb(
	const double	*x,		/* vector */
	const double	*y,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	const double 	beta,		/* scalar */
	double			*z
	);

// B := (A^t)(A) 
void Com_mat_AtA(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	double			*B	/* B[m*m] */
	);

/* z:= Ax */
void Com_mat_Ax(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[m] */
	double			*z		/* return z[n] */
	);

/* z:= A^T x */
void Com_mat_Ax(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[n] */
	double			*z		/* return z[m] */
	);

// compute determinant matrix
double determinant(
	double	*B_,
	int		m
	);

// compute eigenvalues and eigenvectors
int Com_eigen(
	const	double	*A,	// [n*n]
	const	int		n,
	double			*lam,	// [n] :eigenvalues
	double			*V		// [n*n] :eigenvectors
	);

// return true if given values are equal
bool Equal(
	const double	a,
	const double 	b,
	const double	eq
	);



/* ---------------------------------------------------*/
// Template Functions:
/* ---------------------------------------------------*/

// transpose A 
template <typename Type>
void TraMat(
	int 		n,
	int		m,
	Type		*A,	/* [n*m] or [n][m] */
	Type		*T		/* (m,n) */
	);

// print a vector 
template <typename Type>
void printv(
	int	n,
	Type	*x
	);

// print a matrix
template <typename Type>
void printM(
	int	n,
	int	m,
	Type	*x
	);

// generate a zero vector
template <typename Type>
void Gen_ZeroVec(
	int	n,
	Type	*x
	);

// count positive values
template <typename Type>
int	count(
	int	n,		//size
	Type	*x
	);

// bubule sort
template <typename Type>
void Bubble_Sort(
	//input
	const int	n,
	const Type	*x,	// [n]
	//output
	Type			*y,	// [n] or NULL
	int			*z		// [n] or NULL
	);
	

/* ---------------------------------------------------*/
// Inline Functions:
/* ---------------------------------------------------*/

// return ||x||
inline double	Com_nrm(
	const double	*x,		/* vector */
	const int		n			/* size */
	)
{
	int one = 1;
	return cblas_dnrm2( n, x, one);
}


/* copy a vector */
/* y := x */
inline void Copy_vec(
	const double	*x,
	double			*y,
	const int		n		/* size */
	)
{
	int inc=1;	/* increment */
	cblas_dcopy( n, x, inc, y, inc); 
}
// compute a inner product
inline double Com_dot(
	const double	*x,	/* vector */
	const double	*y,	/* vector */
	const int		n		/* size */
	)
{
	int 		one=1;
	return cblas_ddot( n, x, one, y, one);
}

#endif
