#ifndef SCHMIDT_H_
#define SCHMIDT_H_

#include <assert.h>
#include "vector.h"

using namespace std;

class SCHMIDT_M{
	int		n;
	double	min;
	bool		*z;	// [n]

	double	*B_;	// [n*n]
	double	*GS_;	// [n*n]
	double	*nrm;	// [n]	nrm[i] = || (GS_)_i ||^2
	public:

	SCHMIDT_M();											// default constructor	
	SCHMIDT_M( const SCHMIDT_M &source );			// copy constructor
	SCHMIDT_M& operator=( const SCHMIDT_M& );	// assignment operator
	~SCHMIDT_M();										// destructor

	void		setup( const int s_n, const double *s_B_);
	void		compute_GS();

	double	get_min(){
		assert( min > 0.0 );
		return min;
	}

	double	u( const int i, const int j){
		assert( i > j && i < n && j < n);
		assert( nrm[j] >= 0 );
		return (Com_dot( &B_[i*n], &GS_[j*n], n) / nrm[j]);
	}

	void		set_z_i(int	i, bool val){
		assert( z != NULL );
		assert( n > 0 );
		z[i] = val;
	}

	int	get_n(){ return n; }

};

#endif
