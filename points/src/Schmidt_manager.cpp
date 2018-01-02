#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "Schmidt_manager.h"
#include "vector.h"

using namespace std;


#define debug_class	0
#define debug	0

SCHMIDT_M::SCHMIDT_M(){	// default constructor
#if debug_class
	cout << "SCHMIDT_M: default constructor" << endl;
#endif
	n = -1;
	min = -1.0;
	z = NULL;
	B_ = NULL;
	GS_ = NULL;
	nrm = NULL;
}

SCHMIDT_M::SCHMIDT_M( const SCHMIDT_M &source )
{	// copy constructor
#if debug_class
	cout << "SCHMIDT_M: copy constructor" << endl;
#endif
	n = source.n;
	min = source.min;

	if( source.z != NULL ){
		assert( n > 0 );
		z = new bool[n];
		for(int i=0; i<n; i++) z[i] = source.z[i];
	}else{
		z = NULL;
	}

	if( source.B_ != NULL ){
		assert( n > 0 );
		B_ = new double[n*n];
		Copy_vec( source.B_, B_, n*n);
	}else{
		B_ = NULL;
	}

	if( source.GS_ != NULL ){
		assert( n > 0 );
		GS_ = new double[n*n];
		Copy_vec( source.GS_, GS_, n*n);
	}else{
		GS_ = NULL;
	}

	if( source.nrm != NULL ){
		assert( n > 0 );
		nrm = new double[n];
		Copy_vec( source.nrm, nrm, n);
	}else{
		nrm = NULL;
	}

}

// assignment operator
SCHMIDT_M& SCHMIDT_M::operator=( const SCHMIDT_M& source )
{
#if debug_class
	cout << "SCHMIDT_M: assignment operator" << endl;
#endif

	if( this != &source ){ 	
		n = source.n;
		min = source.min;

		if( source.z != NULL ){
			assert( n > 0 );
			z = new bool[n];
			for(int i=0; i<n; i++) z[i] = source.z[i];
		}else{
			z = NULL;
		}
	
		if( source.B_ != NULL ){
			assert( n > 0 );
			B_ = new double[n*n];
			Copy_vec( source.B_, B_, n*n);
		}else{
			B_ = NULL;
		}
	
		if( source.GS_ != NULL ){
			assert( n > 0 );
			GS_ = new double[n*n];
			Copy_vec( source.GS_, GS_, n*n);
		}else{
			GS_ = NULL;
		}
	
		if( source.nrm != NULL ){
			assert( n > 0 );
			nrm = new double[n];
			Copy_vec( source.nrm, nrm, n);
		}else{
			nrm = NULL;
		}
	}
	return *this;
}

// destructor
SCHMIDT_M::~SCHMIDT_M()
{
#if debug_class
	cout << "SCHMIDT_M: destructor" << endl;
#endif
	delete[] z;
	delete[] B_;
	delete[] GS_;
	delete[] nrm;
	z = NULL;
	B_ = NULL;
	GS_ = NULL;
	nrm = NULL;
}

void SCHMIDT_M::compute_GS(){

	assert( GS_ != NULL );
	assert( B_ != NULL );
	assert( n > 0 );
	assert( z != NULL );

	min = -1.0;
	for(int i=0; i<n; i++){
		if( z[i] == true ){
			Copy_vec( &B_[i*n], &GS_[i*n], n);
			for(int j=0; j<i; j++){
				if( z[j] == true ){
					Com_linecomb( &GS_[i*n], &GS_[j*n], n, 1.0, - u(i,j), &GS_[i*n]);
				}
			}
			nrm[i] = Com_dot( &GS_[i*n], &GS_[i*n], n);
			if( min == -1.0 ){
				min = nrm[i];
			}else if( min > nrm[i] ){
				min = nrm[i];
			}
		}else{
			for(int j=0; j<n; j++){
				GS_[(i*n)+j] = 0.0;
			}
			nrm[i] = 0.0;
		}

	}

	//for(int i=0; i<n; i++){
	//	for(int j=0; j<n; j++){
	//		if( i != j) 
	//			cout << Com_dot( &GS_[i*n], &GS_[j*n], n) << endl;
	//	}
	//}
 

}

void SCHMIDT_M::setup(
	const int		s_n,
	const double	*s_B_
	)
{
	assert( s_n > 0 );
	assert( s_B_ != NULL );

	n = s_n;

	B_ = new double[n*n];
	GS_ = new double[n*n];
	nrm = new double[n*n];
	z = new bool[n];

	Copy_vec( s_B_, B_, n*n);

	for(int i=0; i<n; i++){
		z[i] = true;
		nrm[i] = -1.0;
	}
}

