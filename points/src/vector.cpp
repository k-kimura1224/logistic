
#include <iostream>
#include <assert.h>
#include <math.h>

#include "vector.h"

using namespace std;
template	void TraMat( int n, int m, int *A, int *T);
template	void TraMat( int n, int m, double *A, double *T);
template void printv( int n, int *x);
template void printv( int n, double *x);
template void printM( int n, int m, int *x);
template void printM( int n, int m, double *x);
template void Gen_ZeroVec( int n, int *x);
template void Gen_ZeroVec( int n, double *x);
template int count( int n, bool *x);
template void Bubble_Sort( const int n, const double *x, double *y, int *z);

//clapack
extern "C" {
	int dposv_( char* uplo, int* n, int* nrhs, 
	double* A, int* lda, double* x,
	int* ldb, int* info);
}
extern "C" {
	int dgesv_( int *n, int *nrhs, double *a, 
	int *lda, int *ipiv, double *b, int *ldb, int *info);
}
extern "C" {
	int dsysv_( char* uplo, int* n, int* nrhs, 
	double* A, int* lda, int *ipiv, double* b,
	int* ldb, double *work, int *lwork, int* info);
}
extern "C" {
	int dgetrf_( int* m, int* n, double* B_, int* lda, int* ipiv, int* info); 
}
extern "C" {
	int dsyev_( char* jobz, char* uplo, int* N, double* Q, int* lda,
	double* eigenvalues, double* work, int* lwork, int* info);
}

//cblas
#if defined(__APPLE__)
#include <cblas.h>
#else
extern "C" {
	void cblas_dscal( const int N, const double alpha, double *X, const int incX);
}
extern "C" {
	void cblas_daxpy( const int N, const double alpha, const double *x,
	const int incX, double *y, const int incY);
}
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
                         AtlasConj=114};
extern "C" {
void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
}
extern "C" {
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
}
#endif

/* compute Ax=b with Cholesky decomposition */
/* A: positive definite symmetric matrix */
/* if successful, return 0 */
int Com_LS_dposv(
	const double	*A,		/* A[n*n] */
	const double	*b,		/* b[n]	*/
	int		n,			/* size */
	double			*x			/* solution */
	)
{
	char	uplo[2]	=	"L";
	int	one		=	1;
	int	info;
	double	*AA;

	/* alloc */
	AA = new double[n*n];
	
	Copy_vec( A, AA, n*n);
	Copy_vec( b, x, n);

	dposv_( uplo, &n, &one, AA, &n, x, &n, &info);

	/* free */
	delete[] AA;
	
	return info;

}


/* compute Ax=b with LU decomposition */
/* if successful, return 0 */
int Com_LS_dgesv(
	const double	*A,		/* A[n*n] */
	const double	*b,		/* b[n]	*/
	int		n,			/* size */
	double			*x			/* solution */
	)
{
	int	one		=	1;
	int	info;
	double	*AA;
	int		*ipiv;

	/* alloc */
	AA = new double[n*n];
	ipiv = new int[n];
	
	Copy_vec( A, AA, n*n);
	Copy_vec( b, x, n);

	dgesv_( &n, &one, AA, &n, ipiv, x, &n, &info);

	/* free */
	delete[] AA;
	delete[] ipiv;
	
	return info;

}

/* compute Ax=b with Cholesky decomposition */
/* A: symmetric matrix */
/* if successful, return 0 */
int Com_LS_dsysv(
	const double	*A,		/* A[n*n] */
	const double	*b,		/* b[n]	*/
	int		n,			/* size */
	double			*x			/* solution */
	)
{
	char	uplo[2]	=	"L";
	int	one		=	1;
	int	info;
	double	*AA;
	int		*ipiv;
	double	*work;
	int		lwork = n*10;

	/* alloc */
	AA = new double[n*n];
	ipiv = new int[n];
	work = new double[lwork];
	
	Copy_vec( A, AA, n*n);
	Copy_vec( b, x, n);

	dsysv_( uplo, &n, &one, AA, &n, ipiv, x, &n, work, &lwork, &info);

	/* free */
	delete[] AA;
	delete[] ipiv;
	delete[] work;
	
	return info;

}

/* y := (alpha)x */
void Com_scal(
	const double	*x,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	double			*y
	)
{
	int one = 1;
	Copy_vec( x, y, n);

	cblas_dscal( n, alpha, y, one);

}

// compute linear combination
// z := (alpha)x + (beta)y 
void Com_linecomb(
	const double	*x,		/* vector */
	const double	*y,		/* vector */
	const int		n,			/* size */
	const double	alpha,	/* scalar */
	const double 	beta,		/* scalar */
	double			*z
	)
{
	int one=1;

	double	*rslt;
	rslt = new double[n];

	Com_scal( y, n, beta, rslt);

	cblas_daxpy( n, alpha, x, one, rslt, one);

	Copy_vec( rslt, z, n);

	/* free */
	delete[] rslt;

}

/* B := (A^t)(A) */
void Com_mat_AtA(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	double			*B	/* B[m*m] */
	)
{
	double	zero=0.0;
	double	one=1.0;
	

	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
					m, m, n, one, A, n, A, n, zero, B, m);

}

/* z:= Ax */
void Com_mat_Ax(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[m] */
	double			*z		/* return z[n] */
	)
{

	assert( A != NULL );
	assert( n > 0 );
	assert( m > 0 );
	assert( x != NULL );
	assert( z != NULL );

	int one=1;

	if( x == z ){
		double *z_buf = new double[n];
		Gen_ZeroVec( n, z_buf);
		cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, 1.0, A, n, x, one, 0.0, z_buf, one);
		Copy_vec( z_buf, z, n);
		delete[] z_buf;
	}else{
		cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, 1.0, A, n, x, one, 0.0, z, one);
	}
}

/* z:= A^T x */
void Com_mat_Atx(
	const double	*A,	/* A[n*m] */
	const int		n,		/* row */
	const int		m,		/* column */
	const double	*x,	/* x[n] */
	double			*z		/* return z[m] */
	)
{

	assert( A != NULL );
	assert( n > 0 );
	assert( m > 0 );
	assert( x != NULL );
	assert( z != NULL );

	int one=1;

	if( x == z ){
		double *z_buf = new double[m];
		Gen_ZeroVec( m, z_buf);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, x, one, 0.0, z_buf, one);
		Copy_vec( z_buf, z, m);
		delete[] z_buf;
	}else{
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, x, one, 0.0, z, one);
	}
}
// Compute determinant B
double	determinant(
	double			*B_,	// [m*m]
	int				m
	)
{
	int		i;
	double	det;

	// for det
	int	lda = m;
	int	info;
	int	*ipiv = new int[m];
	double *copyB = new double[m*m];

	Copy_vec( B_, copyB, m*m);

	dgetrf_( &m, &m, copyB, &lda, ipiv, &info); 
	
	det=1.0;
	
	for(i=0; i < m-1; i++){ 
		if(ipiv[i] != i+1)  det = -det;
	}
	
	for(i=0; i < m; i++)  det *= copyB[i+m*i];
	
	delete[] ipiv;
	delete[] copyB;

	return  det;
}


// compute eigenvalues and eigenvectors
int Com_eigen(
	const	double	*A,	// [n*n]
	const	int		n,
	double			*lam,	// [n] :eigenvalues
	double			*V		// [n*n] :eigenvectors
	)
{
	assert( A != NULL );
	assert( n > 0 );
	assert( lam != NULL );
	assert( V != NULL );
	assert( A != V );

	Copy_vec( A, V, n*n );

	char	jobz[2]	=	"V";
	char	uplo[2]	=	"U";
	int	N		=	n;
	int	lda	=	n;
	double	*work;
	work = new double[N*3];
	int	lwork	=	N*3;
	int	info;

	dsyev_( jobz, uplo, &N, V, &lda, lam, work, &lwork, &info);

	delete[] work;

	return info;
}

// return true if given values are equal
bool Equal(
	const	double	a,
	const	double	b,
	const	double	ep
	)
{
	bool rslt = false;
	double max = fabs( a );

	if( max < fabs( b ) )	max = fabs( b );
	if( max < 1.0 ) 	max = 1.0;

	if( fabs( a - b ) < max * ep ) rslt = true;

	return rslt;
}

/* ---------------------------------------------------*/
// Template Functions:
/* ---------------------------------------------------*/

// generate a zero vector
template <typename Type>
void Gen_ZeroVec(
	int	n,
	Type	*x
	)
{
	for(int i=0; i<n; i++){
		x[i] = 0;
	}
}

// transpose A 
template <typename Type>
void TraMat(
	int 		n,
	int		m,
	Type		*A,	/* [n*m] or [n][m] */
	Type		*T		/* (m,n) */
){
	int i,j;

	for(i=0; i<m; ++i){
		for(j=0; j<n; ++j){
			*(T+(n*i)+j) = *(A+(j*m)+i);
		}
	}
}

template <typename Type>
void printv(
	int	n,
	Type	*x
	)
{
	for(int i=0; i<n; i++){
		cout << x[i] << " ";
	}
	cout << endl;

}

template <typename Type>
void printM(
	int	n,
	int	m,
	Type	*x
	)
{
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			cout << x[j+(i*m)] << " ";
		}
		cout << endl;
	}

}

template <typename Type>
int count(
	int	n,
	Type	*x
	)
{
	if( x == NULL ) return 0;

	int ct=0;	
	for(int i=0; i<n; i++){
		if( x[i] ) ct++;
	}
	return ct;
}

// bubule sort
template <typename Type>
void Bubble_Sort(
	//input
	const int	n,
	const Type	*x,	// [n]
	//output
	Type			*y,	// [n] or NULL
	int			*z		// [n] or NULL
	)
{
	assert( n > 0 );
	assert( x != NULL );
	assert( y != NULL || z != NULL );

	Type	*a;
	a = new Type[n];	

	// copy
	for(int i=0; i<n; i++){
		a[i] = x[i];
	}

	if( z != NULL ){
		for(int i=0; i<n; i++){
			z[i] = i;
		}
	}

	Type	buf1;
	int	buf2;

	// bubble sort
	for(int i=0; i<(n-1); ++i){
		for(int j=(n-1); j>i; --j){
			if( a[j]>a[j-1] ){
				buf1	=	a[j];
				a[j]	=	a[j-1];
				a[j-1]=	buf1;

				if( z != NULL ){
					buf2	=	z[j];
					z[j]	=	z[j-1];
					z[j-1]=	buf2;
				}
			}
		}
	}

	if( y != NULL ){
		for(int i=0; i<n; i++){
			y[i] = a[i];
		}
	}

	delete[] a;

}
